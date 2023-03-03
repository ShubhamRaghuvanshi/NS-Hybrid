!!          adbsh_serial.f90
!! This subroutine updates the velocity  field in the spectral 
!! code through adam-bashforth iterations. This integrates the
!! linear part exactly and the the contribution from the
!! non-linear part is approximated through adam-bashforth method.
!! This a the most recent modification. This requires one less
!! vector array than the previous one. This is a serial code.
!! This does not use any fourier transform library calls by 
!! itself thus should go unchanged from one serial machine to
!! another.
!! date : April 17th 2003 updated with slaved adbsh scheme.  
!! -----------------------------------------------------------
		subroutine adbsh_serial
	use mod_serial_fluid
	implicit none
	double precision :: temph,rk2,rk2inv,tmr1,tmc1,tmr2,tmc2,tmr3,tmc3, & 
			    p11,p22,p33,p31,p23,p12,maxv,hypvis,alpha,tfact
	integer :: i1,i2,i3,k1,k2,k3,ksqr,ireal,iimag
!! --------------------------------------------------------------
!! step I : add the contribution from the VXW at (t - deltat) 
!! to v(t) and store it in the array containing VXW at (t-deltat)
!! viz the arrays VWk1,VWk2,VWk3.
!! ---------------------------------------------------------------
!	write(*,*) 'be-ab:',VWk1(8,9,12),Vk2(8,9,12),VWk3(8,9,12)
!$omp parallel &
!$omp shared(VWk1,VWk2,VWk3,Vk1,Vk2,Vk3,time_increment,factor,nalias_sqr,vis,vis2,n2,n3,n1hf)&
!$omp private(i1,i2,i3,k1,k2,k3,ksqr,ireal,iimag,hypvis,rk2,alpha,temph,tfact,rk2inv,p11,p22,p33,p12,p23,p31,tmr1,tmr2,tmr3,tmc1,tmc2,tmc3) 
!$omp do

	do i3 = 1,n3
		k3 = (i3-1) - n3*(i3/(n1hf+1))	
		do i2 = 1,n2
			k2 = (i2-1) - n2*(i2/(n1hf+1))	
			do i1 =1,n1hf
				k1 = i1 -1
				ksqr = k1*k1 + k2*k2 + k3*k3 
				ireal = 2*i1-1
				iimag = 2*i1
				hypvis = vis + vis2*rk2
				rk2 = dfloat(ksqr)*factor*factor
				alpha = hypvis*rk2
				temph = time_increment(ksqr)
				tfact=(1.0d0 - temph)/alpha
!! --- -----------DE-ALIAS HERE ------------------------
				if(ksqr.ge.nalias_sqr)then
				VWk1(ireal,i2,i3)= 0.0d0
				VWk1(iimag,i2,i3)= 0.0d0 
!!
				VWk2(ireal,i2,i3)= 0.0d0
				VWk2(iimag,i2,i3)= 0.0d0
!! 
				VWk3(ireal,i2,i3)= 0.0d0
				VWk3(iimag,i2,i3)= 0.0d0
				else
!! ------calculate the projection operator --------------------
				if(ksqr.gt.0) then
				rk2inv=1.0d0/dfloat(ksqr)
				p11=1-k1*k1*rk2inv
				p22=1-k2*k2*rk2inv
				p33=1-k3*k3*rk2inv
				p12= -k1*k2*rk2inv
				p23= -k2*k3*rk2inv
				p31= -k3*k1*rk2inv
				else
				p11=2.0/3.0
				p22=p11
				p33=p22
				p12= -1.0/3.0
				p23=p12
				p31=p23
				tfact = 0.0d0 
				endif
!! -contribution from VXW at (t - deltat) ------
				tmr1 = -0.50d0*VWk1(ireal,i2,i3)*tfact
				tmc1 = -0.50d0*VWk1(iimag,i2,i3)*tfact
!!
				tmr2 = -0.50d0*VWk2(ireal,i2,i3)*tfact
				tmc2 = -0.50d0*VWk2(iimag,i2,i3)*tfact
!!
!!				tmr3 = -0.50d0*VWk3(ireal,i2,i3)*temph*temph
				tmr3 = -0.50d0*VWk3(ireal,i2,i3)*tfact
				tmc3 = -0.50d0*VWk3(iimag,i2,i3)*tfact
!! -make the contribution divergenceless  and add contribution
!! from Vk to it ---------
			VWk1(ireal,i2,i3)=(p11*tmr1+p12*tmr2+p31*tmr3) &
                                  + Vk1(ireal,i2,i3)*temph
			VWk1(iimag,i2,i3)= (p11*tmc1+p12*tmc2+p31*tmc3) &
				   + Vk1(iimag,i2,i3)*temph
!!
				VWk2(ireal,i2,i3)=(p12*tmr1+p22*tmr2+p23*tmr3) &
				   + Vk2(ireal,i2,i3)*temph
				VWk2(iimag,i2,i3)=(p12*tmc1+p22*tmc2+p23*tmc3) &
				   + Vk2(iimag,i2,i3)*temph
!! 
				VWk3(ireal,i2,i3)=(p31*tmr1+p23*tmr2+p33*tmr3) &
				   + Vk3(ireal,i2,i3)*temph
				VWk3(iimag,i2,i3)=(p31*tmc1+p23*tmc2+p33*tmc3) &
				   + Vk3(iimag,i2,i3)*temph
!! --------------------------------
			endif
			enddo
		enddo
	enddo
!$omp end do
!$omp end parallel

!	write(*,*)'ab-mid:',Vk1(5,6,6),VWk1(5,6,6)
!	write(*,*) 'mid-ab:',VWk1(8,9,12),VWk2(8,9,12),VWk3(8,9,12)
!	maxv = maxval(abs(VWk1))
!	write(*,*)'after 1st part',maxv
!! VWk now contains a contribution from both Vk and also VWk 
!! multiplied by the appropriate factors to update Vk ----------
!! -------------------------------------------------------------
!! step II : evaluate VXW in fourier space. To do this first find
!! W by kXV(k). Here we use a new array Wk. Then do in place 
!! inv fft to get Omega in real space. Then do in place inv fft
!! to get V in real space. Then do V X Omega in real space, and 
!! store the result in Omega, i.e. the array Wk. Then do 
!! in place fft to get VXW in fourier space which is again stored
!! in the array Wk. These are done inside the subroutine 
!! nonlin_serial. 
!! --------------------------------------------------------------
	call nonlin_serial
!! --------------------------------------------------------------
!! step III : Now add the contribution of the non-linear part
!! i.e. VXW to the linear part which has been evaluated in 
!! step I and stored in the arrays VWk. Remember at this
!! stage the array Vk contains the real space velocities, and 
!! should be set to zero at this point. Then the contributions
!! from step I and II are to be added to it.If one
!! wants to calculate the real space structure functions then
!! this is the time to calculate it before it is set to zero. 
!! --------------------------------------------------------------
!	call serial_stfunc(Vk1,Vk2,Vk3,n1,S_p,pth_order)
!! --------------------------------------------------------------
	Vk1 = 0.0d0
	Vk2 = 0.0d0
	Vk3 = 0.0d0
!! ----the second stage of updating ----  
!$omp parallel &
!$omp shared(VWk1,VWk2,VWk3,Vk1,Vk2,Vk3,Wk1,Wk2,Wk3,time_increment,factor,nalias_sqr,vis,vis2,n2,n3,n1hf)&
!$omp private(i1,i2,i3,k1,k2,k3,ksqr,ireal,iimag,hypvis,rk2,alpha,temph,tfact,rk2inv,p11,p22,p33,p12,p23,p31,tmr1,tmr2,tmr3,tmc1,tmc2,tmc3) 
!$omp do

	do i3 = 1,n3
		k3 = (i3-1) - n3*(i3/(n1hf+1))	
		do i2 = 1,n2
			k2 = (i2-1) - n2*(i2/(n1hf+1))	
			do i1 =1,n1hf
				k1 = i1 -1
				ksqr = k1*k1 + k2*k2 + k3*k3
				ireal = 2*i1-1
				iimag = 2*i1
				rk2 = dfloat(ksqr)*factor*factor
				temph = time_increment(ksqr)
				hypvis = vis + vis2*rk2
				alpha = hypvis*rk2
				tfact=(1.0d0 - temph)/alpha

!! ------------ DE - ALIAS HERE ------------------
				if(ksqr.ge.nalias_sqr) then
				Wk1(ireal,i2,i3)= 0.0d0
				Wk1(iimag,i2,i3)= 0.0d0
!!		
				Wk2(ireal,i2,i3)= 0.0d0
				Wk2(iimag,i2,i3)= 0.0d0
!! 
				Wk3(ireal,i2,i3)= 0.0d0
				Wk3(iimag,i2,i3)= 0.0d0
				else	
!! ------calculate the projection operator --------------------
				if(ksqr.gt.0) then
				rk2inv=1.0d0/dfloat(ksqr)
				p11=1-k1*k1*rk2inv
				p22=1-k2*k2*rk2inv
				p33=1-k3*k3*rk2inv
				p12= -k1*k2*rk2inv
				p23= -k2*k3*rk2inv
				p31= -k3*k1*rk2inv
				else
				p11=2.0/3.0
				p22=p11
				p33=p22
				p12= -1.0/3.0
				p23=p12
				p31=p23
				tfact=0.0d0 
				endif
!! --------------
				tmr1 = 1.50d0*Wk1(ireal,i2,i3)*tfact
				tmc1 = 1.50d0*Wk1(iimag,i2,i3)*tfact
!!
				tmr2 = 1.50d0*Wk2(ireal,i2,i3)*tfact
				tmc2 = 1.50d0*Wk2(iimag,i2,i3)*tfact
!!
				tmr3 = 1.50d0*Wk3(ireal,i2,i3)*tfact
				tmc3 = 1.50d0*Wk3(iimag,i2,i3)*tfact
!! --------------
				Vk1(ireal,i2,i3)=(p11*tmr1+p12*tmr2+p31*tmr3) &
	                         		+ VWk1(ireal,i2,i3)
				Vk1(iimag,i2,i3)=(p11*tmc1+p12*tmc2+p31*tmc3) &
						+ VWk1(iimag,i2,i3)
!!
				Vk2(ireal,i2,i3)=(p12*tmr1+p22*tmr2+p23*tmr3) &
						+ VWk2(ireal,i2,i3)
				Vk2(iimag,i2,i3)=(p12*tmc1+p22*tmc2+p23*tmc3) &
						+ VWk2(iimag,i2,i3)
!! 
				Vk3(ireal,i2,i3)=(p31*tmr1+p23*tmr2+p33*tmr3) &
						+ VWk3(ireal,i2,i3)
				Vk3(iimag,i2,i3)=(p31*tmc1+p23*tmc2+p33*tmc3) &
						+ VWk3(iimag,i2,i3)
!! ---------------------
			endif
			enddo
		enddo
	enddo
!$omp end do
!$omp end parallel

	call apply_mask_on_Vk

!	write(*,*) 'end-ab:',Wk1(4,5,6),Vk2(8,4,5),VWk2(8,9,12)
!	maxv = maxval(abs(Vk1))
!	write(*,*)'after all adbsh',maxv
!! -store the present Wk in VWk to be used in the next iteration---
	VWk1 = Wk1
	VWk2 = Wk2
	VWk3 = Wk3
!! -----------------------------------------------------
			end subroutine adbsh_serial 
