	!!        global_array.f90
	!! This subroutine evaluates the global arrays used in the serial 
	!! spectral code. This doesnot use any system call thus should
	!! remain unchanged from one serial machine to another. 
	!! changed for slaved adam-bashforth scheme. 
	!! -----------------------------------------------------------------
	subroutine global_array
	use mod_serial_fluid
	implicit none
	!! defination of other variables 
	integer :: l,j,k,k1,k2,k3,ireal,iimag,ksqr,mshl
  integer :: mode_no1,mode_no2,mode_no3,itrial1,jtrial
  integer :: jmin, jmax,  kmin, kmax,  lmin, lmax             
	real*8 	:: rk,rk2,hypvis
	!! --the density of states calculations ----
	!! --the formula below does the folding in from i to j 
	!! where i varies from 1 to n1 and j should be from -n1/2 + 1 
	!! to + n1/2 : 
	!!          j = (i -1) - n1*(i/(n1hf+1))
	!! this has been checked in the program modulo.f in home/res/fortran
	!! in the theory machines.  
	!! -----------------------------------------
	den_state = 0
	time_increment = 0.0d0
	
	do k=n3_low,n3_high
		k3 = (k-1) - n3*(k/(n1hf+1))	
		do j=1,n2
			k2 = (j-1) - n2*(j/(n1hf+1))	
			do l=1,n1hf
				k1=l-1
				!! -------------
				ksqr = k1*k1 + k2*k2 + k3*k3
				rk2 = dfloat(ksqr)*factor*factor
				rk = dsqrt(rk2)
				mshl = nint(rk)
				hypvis = vis + vis2*rk2
				time_increment(ksqr)= dexp(-delta*hypvis*rk2)
				!! -------------
				if((k1.eq.0).or.(k1.eq.n1h)) then
					den_state(mshl)=den_state(mshl) + 1
				else
					den_state(mshl)=den_state(mshl) + 2
				endif
			enddo
		enddo
	enddo
	!! collect the local array values from all the processs (check needed)
	call MPI_Barrier(MPI_COMM_WORLD, ierror) 
	call MPI_Allreduce(den_state , den_state_total, nshell+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierror)
	call MPI_Allreduce(time_increment , time_increment_total, n1d+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierror)
	den_state = den_state_total
	time_increment = time_increment_total
	
	!! ----now list the forced modes  -----
	mode_no1 = den_state(1)
	mode_no2 = den_state(2)	
	mode_no3 = den_state(3)

	allocate(forced_modes1(mode_no1,4),forced_modes2(mode_no2,4), & 
                                    forced_modes3(mode_no3,4))

	!write(*,*) 'GLOBAL', local_n3, size(forced_modes1), mode_no1*4, mode_no2*4, mode_no3*4
	if (ThisTask .eq. 0) then 
		forced_modes1 = 0.0d0
		forced_modes2 = 0.0d0
		forced_modes3 = 0.0d0
		!! -------------------------------------------------------
		count1 = 0 
		count2 = 0
		count3 = 0
	
		kmin = nint((-4.0+1.0)/(1.0-real(n3)/(real(n1hf+1))))
		kmax = nint((4.0+1.0)/(1.0-real(n3)/(real(n1hf+1))))
	
		jmin = nint((-4.0+1.0)/(1.0-real(n2)/(real(n1hf+1))))
		jmax = nint((4.0+1.0)/(1.0-real(n3)/(real(n1hf+1))))

		lmin = -4+1
		lmax = 4+1	
	
	do k=1,n3
		k3 = (k-1) - n3*(k/(n1hf+1))	
		do j=1,n2
			k2 = (j-1) - n2*(j/(n1hf+1))	
			do l=1,n1hf
					k1=l-1
					!! -------------
					ksqr = k1*k1 + k2*k2 + k3*k3
					rk2 = dfloat(ksqr)*factor*factor
					rk = dsqrt(rk2)
					mshl = nint(rk)
					!! ----------------
					if(mshl.eq.1) then
						count1 = count1 + 1
						forced_modes1(count1,1) = l
						forced_modes1(count1,2) = j
						forced_modes1(count1,3) = k
						if((k1.eq.0).or.(k1.eq.n1h))then
							forced_modes1(count1,4) = 1
						else
							forced_modes1(count1,4) = 2
						endif
					else
					endif

					!! ----------------
					if(mshl.eq.2) then
						count2 = count2 + 1
						forced_modes2(count2,1) = l
						forced_modes2(count2,2) = j
						forced_modes2(count2,3) = k
						if((k1.eq.0).or.(k1.eq.n1/2))then
							forced_modes2(count2,4) = 1
						else
							forced_modes2(count2,4) = 2
						endif
					else
					endif
					!! ----------------
					if(mshl.eq.3) then
						count3 = count3 + 1
						forced_modes3(count3,1) = l
						forced_modes3(count3,2) = j
						forced_modes3(count3,3) = k
						if((k1.eq.0).or.(k1.eq.n1/2))then
							forced_modes3(count3,4) = 1
						else
							forced_modes3(count3,4) = 2
						endif
					else
					endif
					!! -------------
				enddo
			enddo
		enddo
	end if 
	call MPI_Barrier(MPI_COMM_WORLD, ierror)

	call MPI_Bcast(count1 , 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)
	call MPI_Bcast(count2 , 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)
	call MPI_Bcast(count3 , 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)

	call MPI_Bcast(forced_modes1 , mode_no1*4, MPI_INT, 0, MPI_COMM_WORLD, ierror)
	call MPI_Bcast(forced_modes2 , mode_no1*4, MPI_INT, 0, MPI_COMM_WORLD, ierror)
	call MPI_Bcast(forced_modes3 , mode_no1*4, MPI_INT, 0, MPI_COMM_WORLD, ierror)


	
	!! -------------checks ----------------------------------------
	if (ThisTask .eq. 0) then 
		open(unit=321,file=trim(cur_dir)//'/density_of_states.out',status='unknown')
		do itrial1 = 0,nshell
			write(321,*)itrial1,den_state(itrial1)
		enddo
		write(321,*)count1,count2,count3
		write(321,*) 'mode 1'
		do itrial1 = 1,count1
			write(321,*) itrial1,(forced_modes1(itrial1,jtrial),jtrial=1,4)
		enddo
		write(321,*) 'mode 2'
		do itrial1 = 1,count2
			write(321,*) itrial1,(forced_modes2(itrial1,jtrial),jtrial=1,4)
		enddo
		write(321,*) 'mode 3'
		do itrial1 = 1,count3
			write(321,*) itrial1,(forced_modes3(itrial1,jtrial),jtrial=1,4)
		enddo
		close(321)	
	end if 
	!! -------------------------------------------------------------------			
	end subroutine global_array
	
	
	
	
