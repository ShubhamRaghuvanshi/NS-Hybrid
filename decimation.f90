subroutine create_mask
  use mod_serial_fluid
	implicit none
  real*8 :: rk2,rk,prob,ran
  integer :: i1,i2,i3,k1,k2,k3,ksqr,mshl,ireal,iimag,i2n,i3n
	
!!  if (nrun .eq. 1) then

	  do i3 = 1,n3
		  k3 = (i3-1) - n3*(i3/(n1hf+1))  
		  do i2 = 1,n2
			  k2 = (i2-1) - n2*(i2/(n1hf+1))  
	      do i1 =1,n1hf
		      k1 = i1 -1
			    ksqr = k1*k1 + k2*k2 + k3*k3 
				  rk2 = factor*factor*dfloat(ksqr)
	        rk = dsqrt(rk2)
		      mshl = nint(rk)
			    ireal = 2*i1 -1
	        iimag = 2*i1
	
					prob = (dfloat(mshl))**(D - 3.0d0)
					call random_number(ran)
					if (ran<prob) then
						mask(ireal,i2,i3) = 1
						mask(iimag,i2,i3) = 1
						N_modes(mshl) = N_modes(mshl) + 1
					else
						mask(ireal,i2,i3) = 0
						mask(iimag,i2,i3) = 0
					endif				

				enddo
			enddo
		enddo				 

	  i1 =1
		k1 = i1 -1
		ireal = 2*i1 -1
		iimag = 2*i1
	  do i3 = 1,n3
		  k3 = (i3-1) - n3*(i3/(n1hf+1))  
		  do i2 = 1,n1hf-1
			  k2 = i2-1   
				if ((k3 > 0) .and. (k3 <n3/2)) then	
					i2n = n2-k2+1			
					i3n = n3-k3+1
					mask(ireal,i2n,i3n) = mask(ireal,i2,i3)
					mask(iimag,i2n,i3n) = mask(iimag,i2,i3)									
				endif
				if (k3 <= 0) then	
					i2n = n2-k2+1			
					i3n = -k3+1
					mask(ireal,i2n,i3n) = mask(ireal,i2,i3)
					mask(iimag,i2n,i3n) = mask(iimag,i2,i3)									
				endif
			enddo
		enddo	
			 
		open(unit = 125,file='den_stat_deci.dat',status='unknown') 
		do i1 = 1,nshell
			write(125,*) i1,N_modes(i1)
		enddo
		close(125)

	  open(unit=122,file='mask.out',form='unformatted',status='unknown')
		write(122)(((mask(i1,i2,i3), i1=1,n1+2),i2=1,n2),i3=1,n3)
	  close(122)

!!	else
!!	  open(unit=122,file='mask.out',form='unformatted',status='unknown')
!!		read(122)(((mask(i1,i2,i3), i1=1,n1+2),i2=1,n2),i3=1,n3)
!!	  close(122)
!!	endif

end subroutine create_mask
!------------------------------------------------

subroutine apply_mask_on_Vk
  use mod_serial_fluid
	implicit none
  integer :: i1,i2,i3,k1,k2,k3,ireal,iimag
	real*8:: tmr1, tmc1, tmr2, tmc2, tmr3, tmc3, mskr, mskc

!$omp parallel &
!$omp shared(Vk1,Vk2,Vk3,mask,n2,n3,n1hf)&
!$omp private(i1,i2,i3,k1,k2,k3,ireal,iimag,tmr1,tmc1,tmr2,tmc2,tmr3,tmc3,mskr,mskc) 
!$omp do	
  do i3 = 1,n3
    k3 = (i3-1) - n3*(i3/(n1hf+1))  
    do i2 = 1,n2
      k2 = (i2-1) - n2*(i2/(n1hf+1))  
      do i1 =1,n1hf
        k1 = i1 -1
        ireal = 2*i1 -1
        iimag = 2*i1

				mskr = dfloat(mask(ireal,i2,i3))				
				mskc = dfloat(mask(iimag,i2,i3))
 				
				tmr1 = Vk1(ireal,i2,i3)
				tmc1 = Vk1(iimag,i2,i3)
				                       
				tmr2 = Vk2(ireal,i2,i3)
				tmc2 = Vk2(iimag,i2,i3)
                               
				tmr3 = Vk3(ireal,i2,i3)
				tmc3 = Vk3(iimag,i2,i3)			
	
				Vk1(ireal,i2,i3) = tmr1*mskr
				Vk1(iimag,i2,i3) = tmc1*mskc
                               
				Vk2(ireal,i2,i3) = tmr2*mskr
				Vk2(iimag,i2,i3) = tmc2*mskc
                               
				Vk3(ireal,i2,i3) = tmr3*mskr
				Vk3(iimag,i2,i3) = tmc3*mskc

			enddo
		enddo
	enddo
!$omp end do
!$omp end parallel

end subroutine apply_mask_on_Vk
!---------------------------------------------------

subroutine apply_mask_on_nonlin
  use mod_serial_fluid
	implicit none
  integer :: i1,i2,i3,k1,k2,k3,ireal,iimag
	real*8:: tmr1, tmc1, tmr2, tmc2, tmr3, tmc3, mskr, mskc
	
!$omp parallel &
!$omp shared(Wk1,Wk2,Wk3,mask,n2,n3,n1hf)&
!$omp private(i1,i2,i3,k1,k2,k3,ireal,iimag,tmr1,tmc1,tmr2,tmc2,tmr3,tmc3,mskr,mskc) 
!$omp do	
  do i3 = 1,n3
    k3 = (i3-1) - n3*(i3/(n1hf+1))  
    do i2 = 1,n2
      k2 = (i2-1) - n2*(i2/(n1hf+1))  
      do i1 =1,n1hf
        k1 = i1 -1
        ireal = 2*i1 -1
        iimag = 2*i1
				
				mskr = dfloat(mask(ireal,i2,i3))				
				mskc = dfloat(mask(iimag,i2,i3))
 				
				tmr1 = Wk1(ireal,i2,i3)
				tmc1 = Wk1(iimag,i2,i3)
				                       
				tmr2 = Wk2(ireal,i2,i3)
				tmc2 = Wk2(iimag,i2,i3)
                               
				tmr3 = Wk3(ireal,i2,i3)
				tmc3 = Wk3(iimag,i2,i3)			
	
				Wk1(ireal,i2,i3) = tmr1*mskr
				Wk1(iimag,i2,i3) = tmc1*mskc
                               
				Wk2(ireal,i2,i3) = tmr2*mskr
				Wk2(iimag,i2,i3) = tmc2*mskc
                               
				Wk3(ireal,i2,i3) = tmr3*mskr
				Wk3(iimag,i2,i3) = tmc3*mskc

			enddo
		enddo
	enddo
!$omp end do
!$omp end parallel

end subroutine apply_mask_on_nonlin
