	!!       serial_rnkt.f90
	!! This is the code for serial runge-kutta iterations for 
	!! the spectral code. This does not make any explicit call to
	!! the fourier transform libraries, thus may go unchanged from
	!! one serial machine to another. 
	!! --------------------------------------------------------------
	subroutine rnkt_serial
	use mod_serial_fluid
	implicit none
	real*8 :: rk2,rk2inv,delta1,p11,p22,p33,p12,p23,p31,tr1,tc1,tr2,&
			tc2,tr3,tc3,veff
	integer :: i1,i2,i3,k1,k2,k3,ksqr,ireal,iimag
	!! ---------------------------------------------------------
	allocate(tm1(n1+2,n2,local_n3),tm2(n1+2,n2,local_n3),tm3(n1+2,n2,local_n3))
	tm1 = 0.0d0
	tm2 = 0.0d0
	tm3 = 0.0d0
	delta1 = delta/2.0d0
	!! ---predictor step ---------------------
	!! first  Vk is stored in the temporary
	!! arrays. tm1,tm2,tm3. 
	!! ------------------------------------
	tm1 = Vk1
	tm2 = Vk2
	tm3 = Vk3
	!! -evaluate the non-linear part ----------------------------
	

	call nonlin_serial


	!! store the VXW in fourier space for first adbsh step 
	VWk1 = Wk1
	VWk2 = Wk2
	VWk3 = Wk3
	!! ------ Vk1 are now real space velocities, but here we bring 
	!! write back the tm1 to Vk1 -------------------------------
	Vk1 = tm1
	Vk2 = tm2
	Vk3 = tm3
	!! -----------------------------------------------------------

	do i3 = n3_low,n3_high
		indx3 = i3 - n3_low + 1
		k3 = (i3-1) - n3*(i3/(n1hf+1))	
		do i2 = 1,n2
			k2 = (i2-1) - n2*(i2/(n1hf+1))	
			do i1 =1,n1hf
				k1 = i1 -1
				ksqr = k1*k1 + k2*k2 + k3*k3 
				rk2 = factor*factor*dfloat(ksqr)
				ireal = 2*i1-1
				iimag=2*i1
!!--------   DE-ALIAS here ---------------------
				if(ksqr.ge.nalias_sqr)then
					Vk1(ireal,i2,indx3)= 0.0d0
		      Vk1(iimag,i2,indx3)= 0.0d0

		      Vk2(ireal,i2,indx3)= 0.0d0
		      Vk2(iimag,i2,indx3)= 0.0d0

		      Vk3(ireal,i2,indx3)= 0.0d0
		      Vk3(iimag,i2,indx3)= 0.0d0
				else
!! -------the projection operator --------------				
					if(rk2.gt.1.0E-4) then
						rk2inv= 1.0d0/dfloat(ksqr)
						p11=1.0d0-k1*k1*rk2inv
						p22=1.0d0-k2*k2*rk2inv
						p33=1.0d0-k3*k3*rk2inv
						p12= -k1*k2*rk2inv
						p23= -k2*k3*rk2inv
						p31= -k3*k1*rk2inv
					else
						p11=2.0d0/3.0d0
						p22=p11
						p33=p22
						p12= -1.0d0/3.0d0
						p23=p12
						p31=p23
					endif
!! ---------------------------------------
				tr1 = Wk1(ireal,i2,indx3)
				tc1 = Wk1(iimag,i2,indx3)
!! 
				tr2 = Wk2(ireal,i2,indx3)
				tc2 = Wk2(iimag,i2,indx3)
!!
				tr3 = Wk3(ireal,i2,indx3)
				tc3 = Wk3(iimag,i2,indx3)
				veff = (vis + vis2*rk2)*rk2
!! ---time stepping  -----------------------------------
!! tm1 etc contains the velocity at t = 0, and Vk1 is also at t = 0
!! now. After this iteration Vk1 is at t = delta/2 = delta1. 
!! This is the predictor step. The contribution to the corrector
!! step comes from the Vk1 of this stage. The contribution to
!! the corrector stage of this Vk1 and also Vk1 at t = 0 (which is
!! already stored in tm1) is stored in tm1. The scheme is done
!! in this roundabout way becuase whenever I call seial-nonlin
!! to evaluate the nonlinear part of the NS equation, the Vk1 it
!! returns is not velocity at fourier space but at real space. 
!! ------------------------------------------------------------
				Vk1(ireal,i2,indx3)=delta1*(p11*tr1+p12*tr2 + &
					p31*tr3 - veff*Vk1(ireal,i2,indx3))+Vk1(ireal,i2,indx3)
				Vk1(iimag,i2,indx3)=delta1*( p11*tc1+ p12*tc2+ &
					p31*tc3 - veff*Vk1(iimag,i2,indx3))+Vk1(iimag,i2,indx3)
!! ----
				tm1(ireal,i2,indx3) = tm1(ireal,i2,indx3) - &
						delta*veff*Vk1(ireal,i2,indx3) 
				tm1(iimag,i2,indx3) = tm1(iimag,i2,indx3) - &
						delta*veff*Vk1(iimag,i2,indx3) 
!! ------------------------
				Vk2(ireal,i2,indx3)=delta1*(p12*tr1+ p22*tr2+ &
					p23*tr3  - veff*Vk2(ireal,i2,indx3))+Vk2(ireal,i2,indx3)
				Vk2(iimag,i2,indx3)=delta1*( p12*tc1+ p22*tc2+ &
					p23*tc3  - veff*Vk2(iimag,i2,indx3))+Vk2(iimag,i2,indx3)
!! ------
				tm2(ireal,i2,indx3) = tm2(ireal,i2,indx3) - &
						delta*veff*Vk2(ireal,i2,indx3) 
				tm2(iimag,i2,indx3) = tm2(iimag,i2,indx3) - &
						delta*veff*Vk2(iimag,i2,indx3) 
!! ---------------------------
				Vk3(ireal,i2,indx3)=delta1*( p31*tr1+ p23*tr2+ &
					p33*tr3  - veff*Vk3(ireal,i2,indx3))+Vk3(ireal,i2,indx3)
				Vk3(iimag,i2,indx3)=delta1*( p31*tc1+ p23*tc2+ &
					p33*tc3  - veff*Vk3(iimag,i2,indx3))+Vk3(iimag,i2,indx3)
!! ----------
				tm3(ireal,i2,indx3) = tm3(ireal,i2,indx3) - &
						delta*veff*Vk3(ireal,i2,indx3) 
				tm3(iimag,i2,indx3) = tm3(iimag,i2,indx3) - &
						delta*veff*Vk3(iimag,i2,indx3) 
!! -------------------------------------------
			endif
			enddo
		enddo
	enddo
!! -the corrector step ----------------------------------
!! -------evaluating the non-linear part ----------------

	!uncomment
	call nonlin_serial

!! --Vk1 at this stage is velocity in real space not fourier space,
!! the contribution from the velocity in fourier space is already 
!! stored in the arrasys tm1. we set Vk1 to zero here.----------
	Vk1 = 0.0d0
	Vk2 = 0.0d0
	Vk3 = 0.0d0
!! --the corrector time-stepping -------------------------------------
	do i3 = n3_low,n3_high
		indx3 = i3 - n3_low + 1
		k3 = (i3-1) - n3*(i3/(n1hf+1))	
		do i2 = 1,n2
			k2 = (i2-1) - n2*(i2/(n1hf+1))	
			do i1 =1,n1hf
				k1 = i1 -1
				ksqr = k1*k1 + k2*k2 + k3*k3 
				rk2 = factor*factor*dfloat(ksqr)
				ireal = 2*i1-1
				iimag = 2*i1
!! ---------- DE-ALIAS here --------------------
				if(ksqr.ge.nalias_sqr)then
					Vk1(ireal,i2,indx3)= 0.0d0
		      Vk1(iimag,i2,indx3)= 0.0d0
					!!
		      Vk2(ireal,i2,indx3)= 0.0d0
		      Vk2(iimag,i2,indx3)= 0.0d0
					!!
		      Vk3(ireal,i2,indx3)= 0.0d0
		      Vk3(iimag,i2,indx3)= 0.0d0
				else

!! -------the projection operator --------------
				if(ksqr.gt.0) then
					rk2inv= 1.0d0/dfloat(ksqr)
					p11=1.0d0-k1*k1*rk2inv
					p22=1.0d0-k2*k2*rk2inv
					p33=1.0d0-k3*k3*rk2inv
					p12= -k1*k2*rk2inv
					p23= -k2*k3*rk2inv
					p31= -k3*k1*rk2inv
				else
					p11=2.0d0/3.0d0
					p22=p11
					p33=p22
					p12= -1.0d0/3.0d0
					p23=p12
					p31=p23
				endif
!! --------------------------------------------------
				tr1=Wk1(ireal,i2,indx3)
				tc1=Wk1(iimag,i2,indx3)
!!
				tr2=Wk2(ireal,i2,indx3)
				tc2=Wk2(iimag,i2,indx3)
!!
				tr3=Wk3(ireal,i2,indx3)
				tc3=Wk3(iimag,i2,indx3)
!! ----------------------
				veff= (vis+vis2*rk2)*rk2
!! ------------------------------------------------
				Vk1(ireal,i2,indx3) =  tm1(ireal,i2,indx3) &
					        + delta*( p11*tr1+ p12*tr2+ p31*tr3 )  
				Vk1(iimag,i2,indx3) = tm1(iimag,i2,indx3) & 
					        + delta*( p11*tc1+ p12*tc2+ p31*tc3 ) 
!! --------------------
				Vk2(ireal,i2,indx3) =  tm2(ireal,i2,indx3) &
					+ delta*( p12*tr1+ p22*tr2+ p23*tr3  )
				Vk2(iimag,i2,indx3) = tm2(iimag,i2,indx3) &
					+ delta*( p12*tc1+ p22*tc2+ p23*tc3 ) 
!! -----------------------
				Vk3(ireal,i2,indx3) =  tm3(ireal,i2,indx3) &
					+ delta*( p31*tr1+ p23*tr2+ p33*tr3  )
				Vk3(iimag,i2,indx3) = tm3(iimag,i2,indx3) &
					+ delta*( p31*tc1+ p23*tc2+ p33*tc3  )
!! ---------------------------------------------------
			endif
			enddo
		enddo
	enddo
!! -----------------------------------------------------------------
	deallocate(tm1,tm2,tm3)
!! ------------------------
			end subroutine rnkt_serial
