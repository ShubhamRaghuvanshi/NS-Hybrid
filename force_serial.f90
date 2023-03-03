
	subroutine force_serial
	use mod_serial_fluid
	implicit none
	real*8 :: multiply,energy1,energy2,energy3, &
		  rescale3,fact1,fact2, &
		  Omega,rk2inv	 
	double precision :: n_cube
	integer :: j,ireal,iimag,i1,i2,i3,k1,k2,k3,ksqr
	real*8 :: p11,p22,p33,p12,p23,p31,tfact
	!! The forcing terms before and after projection
	real*8 :: fxr,fxi,fyr,fyi,fzr,fzi,pfxr,pfxi,pfyr,pfyi,pfzr,pfzi
!! ----------------------------------------------------------------
	n_cube = dfloat(n1*n2*n3)
!! -----------------------
        fk1 = 0.0d0
        fk2 = 0.0d0
        fk3 = 0.0d0

!! --calculate energy in shell 1 -----
	energy1 = 0.0d0 ! initialise the total energy 

!$omp parallel &
!$omp shared(forced_modes1,Vk1,Vk2,Vk3,count1) &
!$omp private(j,ireal,iimag,i2,i3,multiply) &
!$omp reduction(+ :energy1)
!$omp do
	
	do j=1,count1
	   ireal = 2*forced_modes1(j,1)-1
	   iimag = 2*forced_modes1(j,1)
	   i2 = forced_modes1(j,2)
	   i3 = forced_modes1(j,3)
	   multiply = dfloat(forced_modes1(j,4))	
	   energy1 = energy1+(Vk1(ireal,i2,i3)**2+Vk1(iimag,i2,i3)**2 &
	   		  + Vk2(ireal,i2,i3)**2 + Vk2(iimag,i2,i3)**2 &
			  + Vk3(ireal,i2,i3)**2 + Vk3(iimag,i2,i3)**2)*multiply
	enddo
!$omp end do
!$omp end parallel
!! ---
	energy1 = energy1/(n_cube*n_cube)
!! ------------------

!! --calculate energy in shell 2 -----
	energy2 = 0.0d0 ! initialise 

!$omp parallel &
!$omp shared(forced_modes2,Vk1,Vk2,Vk3,count2) &
!$omp private(j,ireal,iimag,i2,i3,multiply) &
!$omp reduction(+ :energy2)
!$omp do

	do j=1,count2
		ireal = 2*forced_modes2(j,1)-1
		iimag = 2*forced_modes2(j,1)
		i2 = forced_modes2(j,2)
		i3 = forced_modes2(j,3)
		multiply = dfloat(forced_modes2(j,4))
	energy2 = energy2+(Vk1(ireal,i2,i3)**2+Vk1(iimag,i2,i3)**2 &
			+ Vk2(ireal,i2,i3)**2 + Vk2(iimag,i2,i3)**2 &
			+ Vk3(ireal,i2,i3)**2 + Vk3(iimag,i2,i3)**2)*multiply
	enddo
!$omp end do
!$omp end parallel

	energy2 = energy2/(n_cube*n_cube)
!! --------------------------------------

!! Start the forcing step
!! ----rescaling factor --------
	   fact1 = fixed_energy1/(2.0d0*energy1)
	   fact2 = fixed_energy2/(2.0d0*energy2)
!! ------rescale the velocities ---------------------------

!$omp parallel &
!$omp shared(forced_modes1,Vk1,Vk2,Vk3,count1,fact1,fk1,fk2,fk3) &
!$omp private(j,i1,ireal,iimag,i2,i3,k1,k2,k3,ksqr,pfxr,pfxi,pfyr,pfyi,pfzr,pfzi) 
!$omp do

	do j=1,count1 !shell 1 
		i1 = forced_modes1(j,1)
		ireal = 2*i1-1
		iimag = 2*i1
		i2 = forced_modes1(j,2)
		i3 = forced_modes1(j,3)  
	!!!! -----------------------
!!
!!  Calculate the forcing
	  pfxr = fact1*Vk1(ireal,i2,i3)
	  pfxi = fact1*Vk1(iimag,i2,i3)
		
	  pfyr = fact1*Vk2(ireal,i2,i3)
	  pfyi = fact1*Vk2(iimag,i2,i3)
		
	  pfzr = fact1*Vk3(ireal,i2,i3)
	  pfzi = fact1*Vk3(iimag,i2,i3)
!!
!---------writing force amplitude to array---------------
								fk1(ireal,i2,i3) = pfxr
								fk1(iimag,i2,i3) = pfxi

                fk2(ireal,i2,i3) = pfyr
                fk2(iimag,i2,i3) = pfyi

                fk3(ireal,i2,i3) = pfzr
                fk3(iimag,i2,i3) = pfzi

	enddo
!$omp end do
!$omp end parallel
!! -------------------------------------------
!$omp parallel &
!$omp shared(forced_modes2,Vk1,Vk2,Vk3,count2,fact2,fk1,fk2,fk3) &
!$omp private(j,i1,ireal,iimag,i2,i3,k1,k2,k3,ksqr,pfxr,pfxi,pfyr,pfyi,pfzr,pfzi) 
!$omp do

	do j=1,count2 !shell 2 
   	i1 = forced_modes2(j,1)
		ireal = 2*i1-1
		iimag = 2*i1
		i2 = forced_modes2(j,2)
		i3 = forced_modes2(j,3)  
!! -----------------------
!!  Calculate the forcing
	  pfxr = fact2*Vk1(ireal,i2,i3)
	  pfxi = fact2*Vk1(iimag,i2,i3)
		
	  pfyr = fact2*Vk2(ireal,i2,i3)
	  pfyi = fact2*Vk2(iimag,i2,i3)
		
	  pfzr = fact2*Vk3(ireal,i2,i3)
	  pfzi = fact2*Vk3(iimag,i2,i3)

!---------writing force amplitude to array---------------

                fk1(ireal,i2,i3) = pfxr
                fk1(iimag,i2,i3) = pfxi

                fk2(ireal,i2,i3) = pfyr
                fk2(iimag,i2,i3) = pfyi

                fk3(ireal,i2,i3) = pfzr
                fk3(iimag,i2,i3) = pfzi	
	enddo
!$omp end do
!$omp end parallel
!! -------------------------------------------------

			end subroutine force_serial
