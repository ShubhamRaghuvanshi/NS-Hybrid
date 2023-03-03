!! 		serial_stfunc.f90 
!! Calculates structure function in real space. 
	subroutine stfunc_serial
	use mod_serial_fluid
	implicit none
	integer :: ix,iy,iz,ilx,ily,ilz,lsqr,ip,modl
	double precision :: onebylsqr,onebyl,deltavx,deltavy,deltavz,& 
						deltav_mod2,deltav_long,deltav_tran,lsftemp, & 
						tsftemp,lx,ly,lz
!	double precision :: check
!! -----------------------------------------------------
	do iz = 1,n3
		do iy = 1,n2
			do ix = 1,n1
				ilx = ix - ix0
				ily = iy - iy0
				ilz = iz - iz0
!! -----------apply periodic boundary condition here ----------------
				ilx = ilx - n1*(ilx/n1hf)
				ily = ily - n2*(ily/n1hf)
				ilz = ilz - n3*(ilz/n1hf)
!				write(*,*) ilx,ily,ilz
!! ---------- ---------------------------------------------------------
				lx = dfloat(ilx)
				ly = dfloat(ily)
				lz = dfloat(ilz)
				lsqr = ilx*ilx + ily*ily + ilz*ilz
				modl = nint(dsqrt(dfloat(lsqr)))
				if(modl.eq.0)then 
				Sp_long(0,1:pmax) = 0.0d0
				Sp_tran(0,1:pmax) = 0.0d0
				else
				onebylsqr = 1.0d0/dfloat(lsqr)
				onebyl = dsqrt(onebylsqr)
				deltavx = Vk1(ix,iy,iz) - Vk1(ix0,iy0,iz0) 
				deltavy = Vk2(ix,iy,iz) - Vk2(ix0,iy0,iz0) 
				deltavz = Vk3(ix,iy,iz) - Vk3(ix0,iy0,iz0)
				deltav_mod2 = deltavx*deltavx +deltavy*deltavy +deltavz*deltavz 
				deltav_long = (deltavx*lx + deltavy*ly + deltavz*lz)*onebyl
				deltav_tran = dsqrt(deltav_mod2 - deltav_long*deltav_long)
!				check = (ilx*ilx + ily*ily + ilz*ilz)*onebylsqr
				lsftemp = 1.0d0
				tsftemp = 1.0d0
				do ip = 1,pmax
					lsftemp = lsftemp*deltav_long
					tsftemp = tsftemp*deltav_tran	
					Sp_long(modl,ip) = Sp_long(modl,ip)+lsftemp/rden_state(modl)
					Sp_tran(modl,ip) = Sp_tran(modl,ip)+tsftemp/rden_state(modl)
!					write(*,*) check,deltav_tran,modl
				enddo
				endif
!!-------------------------------------------------------------------
			enddo
		enddo
	enddo
!! ------------------------------------------------------------------------
	end subroutine stfunc_serial
					 
 
