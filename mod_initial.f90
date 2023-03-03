	module mod_initial
	use mod_serial_fluid
	use omplib
!	use omplib
	implicit none

	contains

	subroutine read_input_params
		do itask = 0, NTask -1
			call MPI_Barrier(MPI_COMM_WORLD, ierror)
			if (itask .eq. ThisTask ) then 
	  		open(unit=10,file='flu.in',status='old')
  			read(10,*)
  			read(10,*)nn,delta,vis,vis2
				read(10,*)
				read(10,*)maxiter,nrun,navg
				read(10,*)
				read(10,*)forcing,kini,amp,fe3
				read(10,*)
				read(10,*)particl,nprtcl,nstokes,nrunpart
				close(10)
  		end if 
  	end do 	
	endsubroutine read_input_params
!--------------------------------------------------------------------------------------------

	subroutine initialise_variable
  	n1 = nn 
  	n2 = nn
  	n3 = nn
		
		if(mod(n1, NTask) .ne. 0 .or. mod(n2*n3, NTask) .ne. 0  ) then
			write(*,*) "Data dimensions not suitable for FFT, Aborting"
			call endrun(101)	
		end if 

		!! n1d : maximum value k^2 can take
		!! nshell : radius of nth shell = sqrt(k^2) = k =sqrt(n1d)
  	n1d = 3*(nn/2)**2
		!!
  	n1h = n1/2
  	n1hf = n1/2 + 1

		!! ----------
  	alloc_local = fftw_mpi_local_size_3d(n3, n2, n1hf, MPI_COMM_WORLD, &
  	& local_n3, local_n3_offset)

		n3_low  = local_n3_offset + 1
		n3_high	=	local_n3_offset  + local_n3 
		cdatar  = fftw_alloc_real(2 * alloc_local)
		cdatac  = fftw_alloc_complex( alloc_local)	

		!write(*,*) "MOD_INITIAL :", n3_low, n3_high

		call c_f_pointer(cdatar, r_data, [n1+2,n2,local_n3])
		call c_f_pointer(cdatac, c_data, [n1hf,n2,local_n3]) 	
		!! ---------
		
		nshell = int(1.732*nn/2.0d0) + 1 
		nalias_sqr = nn*nn/9
		zi = dcmplx(0.0d0,1.0d0)
		zone = dcmplx(1.0d0,0.0d0)
		!! The length of the box is fixed to be twice pi.  
		pi = 4.0d0*datan(1.0d0)
		length = 2*pi
		factor = 2*pi/length  
		scale = 1.0d0/(dfloat(n1)*dfloat(n2)*dfloat(n3))

		dx=length/dble(n1)
		dy=length/dble(n2)
		dz=length/dble(n3)
		dt_by_dx=delta/dx   

		!! %allocate and calculate the global arrays --------
		allocate(den_state(0:nshell))
		allocate(den_state_total(0:nshell))
		
		allocate(time_increment(0:n1d))
		allocate(time_increment_total(0:n1d))
		
		allocate(N_modes(0:nshell))
				
		den_state = 0
		time_increment = 0.0d0
		N_modes = 0
		!! %---and also energy and dissipation arrays ------
		allocate(E_Omega(0:nshell,2),iniEk(0:nshell))
		allocate(E_Omega_total(0:nshell,2) )
	
		E_Omega = 0.0d0

		!! %-velocity array are allocated  ----
		allocate(Vk1(n1+2,n2,local_n3),Vk2(n1+2,n2,local_n3),Vk3(n1+2,n2,local_n3))
		Vk1 = 0.0d0
		Vk2 = 0.0d0
		Vk3 = 0.0d0
		!! %-allocate other system size arrays too ..
		allocate(VWk1(n1+2,n2,local_n3),VWk2(n1+2,n2,local_n3),VWk3(n1+2,n2,local_n3))
		VWk1 = 0.0d0
		VWk2 = 0.0d0
		VWk3 = 0.0d0
		allocate(Wk1(n1+2,n2,local_n3),Wk2(n1+2,n2,local_n3),Wk3(n1+2,n2,local_n3))
		Wk1 = 0.0d0
		Wk2 = 0.0d0
		Wk3 = 0.0d0

		allocate(mask(n1+2,n2,local_n3))
		mask = 0

		allocate(dx_ux(n1+2,n2,local_n3),dx_uy(n1+2,n2,local_n3),dx_uz(n1+2,n2,local_n3))
		allocate(dy_ux(n1+2,n2,local_n3),dy_uy(n1+2,n2,local_n3),dy_uz(n1+2,n2,local_n3))
		allocate(dz_ux(n1+2,n2,local_n3),dz_uy(n1+2,n2,local_n3),dz_uz(n1+2,n2,local_n3))
		dx_ux = 0.0d0; dx_uy = 0.0d0; dx_uz = 0.0d0 
		dy_ux = 0.0d0; dy_uy = 0.0d0; dy_uz = 0.0d0 
		dz_ux = 0.0d0; dz_uy = 0.0d0; dz_uz = 0.0d0 

		allocate(fk1(n1+2,n2,local_n3),fk2(n1+2,n2,local_n3),fk3(n1+2,n2,local_n3))
		fk1 = 0.0d0
		fk2 = 0.0d0
		fk3 = 0.0d0

		!allocate(uxe(64*64*64),uye(64*64*64),uze(64*64*64))
		!allocate(dx_uxe(n1*n2*local_n3),dx_uye(n1*n2*local_n3),dx_uze(n1*n2*local_n3))
		!allocate(dy_uxe(n1*n2*local_n3),dy_uye(n1*n2*local_n3),dy_uze(n1*n2*local_n3))
		!allocate(dz_uxe(n1*n2*local_n3),dz_uye(n1*n2*local_n3),dz_uze(n1*n2*local_n3))
		!uxe = 0.0d0; uye = 0.0d0; uze = 0.0d0
		!dx_uxe = 0.0d0; dx_uye = 0.0d0; dx_uze = 0.0d0
		!dy_uxe = 0.0d0; dy_uye = 0.0d0; dy_uze = 0.0d0
		!dz_uxe = 0.0d0; dz_uye = 0.0d0; dz_uze = 0.0d0
	
		!! %          peculiar for FFTW                                          *
		!! % We set up the plan for FFTW here. At present we are not using WISDOM*
		dim(1) = n1
		dim(2) = n2
		dim(3) = n3
		!! %-------create plan for forward transform -------------------------

		call dfftw_plan_with_nthreads(nthreads)

		plan_r2c = fftw_mpi_plan_dft_r2c_3d(n3, n2, n1 , r_data, c_data, &
		&         MPI_COMM_WORLD, FFTW_MEASURE)

		plan_c2r = fftw_mpi_plan_dft_c2r_3d(n3, n2, n1, c_data, r_data, &
		&         MPI_COMM_WORLD, FFTW_MEASURE)

		!! %-------------------******************------------------------------
	end subroutine initialise_variable
!------------------------------------------------------------------------------------------

	subroutine open_output_files
		if (ThisTask .eq. 0) then 
			open(unit=200,file='energy.dat',status='unknown')
			open(unit=201,file='params.dat',status='unknown')
			open(unit=202,file='energy_injection.dat',status='unknown')

			!!  open(unit=15,file='euler/uxe.out',status='unknown')
			!!  open(unit=16,file='euler/uye.out',status='unknown')
			!!  open(unit=17,file='euler/uze.out',status='unknown')
		
			!open(unit=18,file='euler/dx_uxe.out',status='unknown')
			!open(unit=19,file='euler/dx_uye.out',status='unknown')
			!open(unit=20,file='euler/dx_uze.out',status='unknown')
			!open(unit=21,file='euler/dy_uxe.out',status='unknown')
			!open(unit=22,file='euler/dy_uye.out',status='unknown')
			!open(unit=23,file='euler/dy_uze.out',status='unknown')
			!open(unit=24,file='euler/dz_uxe.out',status='unknown')
			!open(unit=25,file='euler/dz_uye.out',status='unknown')
			!open(unit=26,file='euler/dz_uze.out',status='unknown')
		end if 			
	end subroutine open_output_files
!--------------------------------------------------------------------------------------------

	subroutine initial_configuration
	implicit none
	real*8,dimension(3) :: ran
	real*8 :: rk2,rk,ek1,vk,p11,p12,p31,p22,p33,p23,rk2inv
	real*8 :: ampli,Vk1_real,Vk1_imag,Vk2_real,Vk2_imag,Vk3_real,Vk3_imag,n_cube
	integer :: i1,i2,i3,k1,k2,k3,ksqr,mshl,ik,iseed,ireal,iimag
	!! -----------------------------------------------------------------
	if (nrun .eq. 1) then	
		ampli=3.0
		iseed =80629.0 + real(3.0*ThisTask)
		n_cube = dfloat(n1)*dfloat(n2)*dfloat(n3)
		indx3 = 0 
    call random_number(ran)
		do i3 = n3_low,n3_high
			indx3 = indx3 + 1
			k3 = (i3-1) - n3*(i3/(n1hf+1))	
			do i2 = 1,n2
				k2 = (i2-1) - n2*(i2/(n1hf+1))	
				do i1 =1,n1hf
					k1 = i1 -1
					ksqr = k1*k1 + k2*k2 + k3*k3 
					rk2 = factor*factor*dfloat(ksqr)
					rk = dsqrt(rk2)
					mshl = nint(rk)
          ek1 = amp*1.0d0*rk2*rk2*exp(-2.0d0*rk2/(kini*kini))
					vk =  dsqrt(ek1/(3*den_state(mshl)))
					vk = vk * n_cube 
					ireal=2*i1-1
					iimag=2*i1
          if(k1.ne.0)then
            call random_number(ran)
          else
            ran = 0.0d0
          endif
					Vk1_real = vk*dcos(2.0d0*pi*ran(1))
					Vk1_imag = vk*dsin(2.0d0*pi*ran(1))
					!! ----------------	
					Vk2_real = vk*dcos(2.0d0*pi*ran(2))
					Vk2_imag = vk*dsin(2.0d0*pi*ran(2))	
					!! ----------------		
					Vk3_real = vk*dcos(2.0d0*pi*ran(3))
					Vk3_imag = vk*dsin(2.0d0*pi*ran(3))
					!! ----------------

					if(ksqr.gt.1.0E-4) then
						rk2inv=1.0d0/dfloat(ksqr)
						p11 = 1.0d0-dfloat(k1*k1)*rk2inv
						p22 = 1.0d0-dfloat(k2*k2)*rk2inv
						p33 = 1.0d0	-dfloat(k3*k3)*rk2inv
						p12 = -dfloat(k1*k2)*rk2inv
						p23 = -dfloat(k2*k3)*rk2inv
						p31 = -dfloat(k3*k1)*rk2inv
					else
						p11=2.0d0/3.0d0
						p22=p11
						p33=p22
						p12= -1.0d0/3.0d0
						p23=p12
						p31=p23 
					endif
					!!-----------------
					Vk1(ireal,i2,indx3)=p11*Vk1_real+p12*Vk2_real+p31*Vk3_real
					Vk1(iimag,i2,indx3)=p11*Vk1_imag+p12*Vk2_imag+p31*Vk3_imag

					Vk2(ireal,i2,indx3)=p12*Vk1_real+p22*Vk2_real+p23*Vk3_real
					Vk2(iimag,i2,indx3)=p12*Vk1_imag+p22*Vk2_imag+p23*Vk3_imag

					Vk3(ireal,i2,indx3)=p31*Vk1_real+p23*Vk2_real+p33*Vk3_real
					Vk3(iimag,i2,indx3)=p31*Vk1_imag+p23*Vk2_imag+p33*Vk3_imag
				
				enddo
			enddo
		enddo

	
		!		call apply_mask_on_Vk
	  call rnkt_serial
	  
	
		!! first go to the CM frame of the fluid 
		if (ThisTask .eq. 0) then 
			Vk1(1,1,1) = 0.0d0
			Vk2(1,1,1) = 0.0d0
			Vk3(1,1,1) = 0.0d0 
			!this sets the real part to be zero,
			!!  the imaginary part should be zero anyway, otherwise the real 
			!! space velocities wont be real. But we still set them to zero.
			Vk1(2,1,1) = 0.0d0
			Vk2(2,1,1) = 0.0d0
			Vk3(2,1,1) = 0.0d0
		end if
	
	else
		do itask = 0, NTask-1
			call MPI_Barrier(MPI_COMM_WORLD, ierror)
			if (itask .eq. ThisTask ) then 
				open(unit=11,file='Vk.in',form='unformatted',status='old')
	  		read(11)(((Vk1(i1,i2,i3),Vk2(i1,i2,i3),Vk3(i1,i2,i3), &
														  i1=1,n1+2),i2=1,n2),i3=1,local_n3)
				close(11)
	  		open(unit=12,file='VXW.in',form='unformatted',status='old')
				read(12)(((VWk1(i1,i2,i3),VWk2(i1,i2,i3),VWk3(i1,i2,i3), &
													     i1=1,n1+2),i2=1,n2),i3=1,local_n3)
				close(12)
			end if 
		end do
	endif

!! ------------------------------------------------------------
	end subroutine initial_configuration

	subroutine fftw_r2c_3D(r_in)
	real*8, dimension(n1+2, n2,  local_n3) :: r_in
	integer*8 :: ix , jx, kx, ir, im 
	
	do kx = 1, local_n3
		do jx = 1, n2
			do ix = 1, n1
				r_data(ix, jx, kx) = r_in(ix, jx, kx)
			end do
		end do 		
	end do

	call fftw_mpi_execute_dft_r2c(plan_r2c, r_data, c_data)

	do kx = 1, local_n3
		do jx = 1, n2
			do ix = 1, n1hf
					ir=2*ix-1
					im=2*ix
					r_in(ir, jx, kx) = real(c_data(ix, jx, kx))
					r_in(im, jx, kx) = aimag(c_data(ix, jx, kx))
			end do
		end do	
	end do
	end subroutine fftw_r2c_3D

!---------------------------------------------------------	

 	subroutine fftw_c2r_3D( r_in)

	real*8, dimension(n1+2, n2,  local_n3) :: r_in 
	!complex(kind=8), dimension(global_N1h, global_N2, local_N3) :: c_in
	integer*8 :: ix , jx, kx, ir, im
		
	do kx = 1, local_n3
		do jx = 1, n2
			do ix = 1, n1hf
				ir=2*ix-1
				im=2*ix
				c_data(ix, jx, kx) = dcmplx ( r_in(ir, jx, kx), r_in(im, jx, kx) )
			end do
		end do 		
	end do

	call fftw_mpi_execute_dft_c2r(plan_c2r, c_data, r_data)

	do kx = 1, local_n3
		do jx = 1, n2
			do ix = 1, n1
				r_in(ix, jx, kx) = r_data(ix, jx, kx)
			end do
		end do 	
	end do
	end subroutine fftw_c2r_3D

!  This function aborts the simulations. If a single processors wants an
!  immediate termination, the function ought to be called with ierr>0. A
!  bunch of MPI-error messages may also appear in this case.  For ierr=0,
!  MPI is gracefully cleaned up, but this requires that all processors
!  call endrun(). All processors may also call it with ierr>0.

subroutine endrun(ierr)
	integer :: ierr
	if(ierr .gt. 0) then 
		call MPI_Abort(MPI_COMM_WORLD, ierr, ierror)
	end if 
	call MPI_Finalize(ierror)
	call exit(0)
end subroutine endrun
	
end module mod_initial






















