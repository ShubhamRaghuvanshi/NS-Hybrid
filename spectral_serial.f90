!! This is a code for Direct Numerical Simulation (DNS) of Navier-Stokes
!! equation, in a periodic box. We simulate a divergence-less flow. 
!! This is a spectral code, i.e. the derivative are taken in Fourier Space. 
!!This is a 2/3rd de-aliased code ,this is used to prevent the blow up
!!of the tail b'cz of not resolving the largeK(or small l) modes(length scales).!!The other technique to avoid this is to add a viscousity(vis2*k^2) term to 
!! the viscousity.
!! ---------------------------------------------------------------------- 
	program spectral_serial
	use mod_serial_fluid
	use mod_initial
	use	omplib
!	use mod_particle
	implicit none
	!! other defination of variables 
	integer :: ispectra,i1,i2,i3,iouter,inner,iter,p
	integer :: itrial1,icall
  real*8 :: En,Omega,sumEkbyk,realk,time !,OMP_get_wtime
  character(1000)::fnum
	character*500::cnpar,formp


	!! % Initialise MPI and fftw_mpi--------------------------------
	call	MPI_INIT(ierror)
	call	MPI_Comm_rank(MPI_COMM_WORLD, ThisTask, ierror)
	call	MPI_Comm_size(MPI_COMM_WORLD, NTask, ierror)
	call  fftw_mpi_init()
	call  dfftw_init_threads(iret)

	!! %---read initial params ----
	call read_input_params

	!---------------------
	if (ThisTask .eq. 0) then 
		!$OMP PARALLEL 
		!$OMP MASTER
 		call get_total_threads(nthreads) 
    
 		!$OMP END MASTER
		!$OMP END PARALLEL
		call system('rm -rf vel spectra euler')
		call system('mkdir vel spectra euler')
		call system('getconf _NPROCESSORS_ONLN  > NPROC')
		open(unit=101, file='NPROC', action='read')
		read(101, *) n_proc
		close(101)
!		write(*,*) 'Number of processes : ',NTask
!		write(*,*) 'Number of threads per process= ',nthreads		
!		if(nthreads*NTask .gt. n_proc) then 
!			write(*,*) 'Maximum number of processes exceeded for this machine, Aborting'
!			call endrun(666) 
!		end if 
	end if

!	call open_output_files


	call initialise_variable


	call global_array


	!	if (particl) call allocate_particle
	!	call create_mask

	call initial_configuration
	
	call energy_serial 

	En = 0.0d0
	Omega = 0.0d0
	if (ThisTask .eq. 0) then 	
		open(unit=110,file='initial_spectra.out',status='unknown')
			do ispectra = 0,nshell
			!write(*,*) ispectra,E_Omega(ispectra,1),E_Omega(ispectra,2)
			write(110,*) ispectra,E_Omega(ispectra,1),E_Omega(ispectra,2)
		enddo
		close(110)
	end if
	En = sum(E_Omega(:,1))
	Omega = sum(E_Omega(:,2))
!---------------
	
	

!	if (particl) then 
!		call initialize_particle
!		call open_output_particle
!	endif

  vel_writer = 0
	iter = 0; cnt = 0      
!	t1=OMP_get_wtime()
	call cpu_time(t1)
	do iouter = 1,maxiter/navg
	  call MPI_Barrier(MPI_COMM_WORLD, ierror)
		count_out=iouter
		do inner = 1,navg
	    iter = iter + 1
      cnt = iter
       	
      if (forcing) call force_serial
			
      call rnkt_serial

!! first go to the CM frame of the fluid 
if (ThisTask .eq. 0) then 
			Vk1(1,1,1) = 0.0d0
			Vk2(1,1,1) = 0.0d0
			Vk3(1,1,1) = 0.0d0 
!this sets the real part to be zero, the imaginary part
!should be zero anyway, otherwise  the real space 
!velocities wont be real. But we still set them to zero.
			Vk1(2,1,1) = 0.0d0
			Vk2(2,1,1) = 0.0d0
			Vk3(2,1,1) = 0.0d0
end if
!! and after each call the forcing function is enforced. 
			
			E_Omega = 0.0d0
			call energy_serial
if (ThisTask .eq. 0) then 			
			En = 0.0d0
			Omega = 0.0d0
			SumEkbyk=0.0d0
			En = sum(E_Omega(:,1))
			Omega = sum(E_Omega(:,2))
			tlr_micro_scale = dsqrt(En/Omega)
			vkrms = dsqrt(2.0d0*En/3.0d0)
!! The Taylor Micro Scale Reynolds number 
			Rlambda = vkrms*tlr_micro_scale/vis
      c_cfl=umax*dt_by_dx                         
			time = (iter*delta*vkrms)/(2.0d0*pi)	
			write(200,*)time,En,Omega  
      write(201,*)Rlambda,c_cfl
			write(*,*)time,En,Omega  
end if 			
		enddo
!! writing the spectra after every navg number of steps
!!------------------------------------------------------
	if (ThisTask .eq. 0) then
		if(mod(cnt,navg).eq.0) then
			write(fnum,'(g8.0)') iouter
			open(unit=224,file='spectra/spectra'//trim(adjustl(fnum))//'.out',status='unknown')
			do ispectra= 0,nshell
				write(224,*)ispectra,E_Omega(ispectra,1),E_Omega(ispectra,2)
			end do
			close(224)
		end if
	end if     
!!------------------------------------------------------
!	open(unit=11,file='vel/Vk'//trim(adjustl(fnum))//'.in',form='unformatted',status='unknown')
!	write(11)(((Vk1(i1,i2,i3),Vk2(i1,i2,i3),Vk3(i1,i2,i3), &
!         i1=1,n1+2),i2=1,n2),i3=1,n3)
!	close(11)
	enddo
	!t2=OMP_get_wtime()
	call cpu_time(t2)
!! First write down the velocity arrays. 
do itask = 0, NTask-1
	call MPI_Barrier(MPI_COMM_WORLD, ierror)
	if (itask .eq. ThisTask ) then      
		open(unit=11,file='Vk.in',form='unformatted',status='unknown')
		write(11)(((Vk1(i1,i2,i3),Vk2(i1,i2,i3),Vk3(i1,i2,i3), &
         i1=1,n1+2),i2=1,n2),i3=1,local_n3)
		close(11)
		!!  and box(t - delta t) from file. 
		open(unit=12,file='VXW.in',form='unformatted',status='unknown')
		write(12)(((VWk1(i1,i2,i3),VWk2(i1,i2,i3),VWk3(i1,i2,i3), &
         i1=1,n1+2),i2=1,n2),i3=1,local_n3)
		close(12)
	end if 
end do 

if (ThisTask .eq. 0) then
	En = 0.0d0
	Omega = 0.0d0
	sumEkbyk = 0.0d0
	write(110,*)
	do ispectra = 1,nshell
		En = En + E_Omega(ispectra,1)
		Omega = Omega + E_Omega(ispectra,2)
		sumEkbyk = sumEkbyk + E_Omega(ispectra,1)/dfloat(ispectra)  
		write(110,*) ispectra,E_Omega(ispectra,1),E_Omega(ispectra,2)
	enddo
	close(110)
!! The total energy and dissipation is then used to calculate 
!! the other averaged quantities which characterise the simulation.
!! The energy dissipation rate $\epsilon = \nu \Omega $         
	edr = vis*Omega
!! The Taylor Micro scale 
	tlr_micro_scale = dsqrt(5.0d0*En/Omega)
!! The dissipation wave number, 
	dissipation_scale = (edr/(vis**3))**(0.25)
!! The root mean square velocity, 
	vkrms = dsqrt(2.0d0*En/3.0d0)
!! The Taylor Micro Scale Reynolds number 
	Rlambda = vkrms*tlr_micro_scale/vis
!! the integral scale, 
	integral_scale = (3.0d0*pi*sumEkbyk)/(4.0d0*En)
!! and the large eddy turnover time. 
	tau_eddy = (integral_scale/vkrms)/delta
!! Defination of Integral Scale is taken from Ref.(\cite{pope}) page 240, 
!! Eq (6.259) and Eq (6.260). 
!! The parameter characterising the simulation are written down in the
!! file parameters.out 
	open(unit=112,file='parameters.out',status='unknown')
	write(112,*) 'viscosity (nu)    :',vis,' hyperviscosity (\nu_h)  :',vis2
	write(112,*) 'energy diss-rate(epsilon)  :',edr
	write(112,*) 'Mean Energy                 :',En
	write(112,*) 'Mean Enstrophy(Omega)      :',Omega
	write(112,*) 'Box-Length                  :',length	
	write(112,*) 'Taylor-Micro-Scale (lambda):',tlr_micro_scale
	write(112,*) 'Dissipation Scale(k_d)      :',dissipation_scale
	write(112,*) 'rms-velocity(u_{rms})       :',vkrms
	write(112,*) 'Re_Lambda (Re_{lambda})    :',Rlambda
	write(112,*) 'Integral Scale(L)           :',integral_scale
	write(112,*) 'large eddy turnover time (tau_{eddy}) :',tau_eddy
	write(112,*) 'no. of itn for tau_eddy',tau_eddy/delta
	write(112,*) 'averaging done over',tau_eddy/(delta*maxiter),'tau_{eddy} times'
	close(112)
	!close(200)
	close(201)
	close(202)
	!close(13)
	close(14)
end if	
!! Finally we deallocate the arrays allocated,  
	deallocate(den_state)
	deallocate(den_state_total)
	deallocate(time_increment)
	deallocate(time_increment_total)
	deallocate(Vk1,Vk2,Vk3)
	deallocate(VWk1,VWk2,VWk3)
	deallocate(Wk1,Wk2,Wk3)
	deallocate(E_Omega)
	deallocate(mask)
	deallocate(dx_ux,dx_uy,dx_uz)
	deallocate(dy_ux,dy_uy,dy_uz)
	deallocate(dz_ux,dz_uy,dz_uz)
	deallocate(fk1,fk2,fk3)
	!deallocate(uxe,uye,uze)
	!deallocate(dx_uxe,dx_uye,dx_uze)
	!deallocate(dy_uxe,dy_uye,dy_uze)
	!deallocate(dz_uxe,dz_uye,dz_uze)

	
!! Destroy the plans created (for FFTW only ) 
	call fftw_destroy_plan(plan_r2c)
	call fftw_destroy_plan(plan_c2r)

!	call rfftwnd_f77_destroy_plan(plan_r2c)
!	call rfftwnd_f77_destroy_plan(plan_c2r)

	t=t2-t1
	call MPI_Allreduce(t , t_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierror)
	if (ThisTask .eq. 0) then 
		print*,"Time taken=",t_all/real(nthreads*NTask)
	end if  
!! %-------------*********************-------------------------------
!! and end our program 		
	call MPI_Finalize(ierror)	
	end program spectral_serial 
	
	
	
	
	
	
