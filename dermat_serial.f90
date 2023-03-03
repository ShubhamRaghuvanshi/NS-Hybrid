        
       
	subroutine dermat_serial
	use mod_serial_fluid 
	use mod_initial
	implicit none
	integer::i1,i2,i3,k1,k2,k3,ireal,iimag
	real*8::Q,R,div
	real*8,dimension(3,3)::der_mat,der_mat2,der_mat3
               
!$omp parallel &
!$omp shared(Vk1,Vk2,Vk3,dx_ux,dx_uy,dx_uz,dy_ux,dy_uy,dy_uz,dz_ux,dz_uy,dz_uz,n2,n3,n1hf)&
!$omp private(i1,i2,i3,indx3,k1,k2,k3,ireal,iimag) 
!$omp do

	do i3 = n3_low,n3_high
		indx3 = i3 - n3_low + 1
		k3 = (i3-1) - n3*(i3/(n1hf+1))  
			do i2 = 1,n2
			k2 = (i2-1) - n2*(i2/(n1hf+1))  
			do i1 =1,n1hf
				k1 = i1 -1
				ireal = 2*i1 -1
				iimag = 2*i1
				!-------------------------------------------------------------------------------
				dx_ux(ireal,i2,indx3)= -k1*Vk1(iimag,i2,indx3)
				dx_ux(iimag,i2,indx3)= k1*Vk1(ireal,i2,indx3)
				!-------------------------------------------------------------------------------
				dy_ux(ireal,i2,indx3)= -k2*Vk1(iimag,i2,indx3)
				dy_ux(iimag,i2,indx3)= k2*Vk1(ireal,i2,indx3)
				!-------------------------------------------------------------------------------
				dz_ux(ireal,i2,indx3)= -k3*Vk1(iimag,i2,indx3)
				dz_ux(iimag,i2,indx3)= k3*Vk1(ireal,i2,indx3)
				!-------------------------------------------------------------------------------
				dx_uy(ireal,i2,indx3)= -k1*Vk2(iimag,i2,indx3)
				dx_uy(iimag,i2,indx3)= k1*Vk2(ireal,i2,indx3)
				!-------------------------------------------------------------------------------
				dy_uy(ireal,i2,indx3)= -k2*Vk2(iimag,i2,indx3)
				dy_uy(iimag,i2,indx3)= k2*Vk2(ireal,i2,indx3)
				!-------------------------------------------------------------------------------
				dz_uy(ireal,i2,indx3)= -k3*Vk2(iimag,i2,indx3)
				dz_uy(iimag,i2,indx3)= k3*Vk2(ireal,i2,indx3)
				!-------------------------------------------------------------------------------
				dx_uz(ireal,i2,indx3)= -k1*Vk3(iimag,i2,indx3)
				dx_uz(iimag,i2,indx3)= k1*Vk3(ireal,i2,indx3)
				!-------------------------------------------------------------------------------
				dy_uz(ireal,i2,indx3)= -k2*Vk3(iimag,i2,indx3)
				dy_uz(iimag,i2,indx3)= k2*Vk3(ireal,i2,indx3)
				!-------------------------------------------------------------------------------
				dz_uz(ireal,i2,indx3)= -k3*Vk3(iimag,i2,indx3)
				dz_uz(iimag,i2,indx3)= k3*Vk3(ireal,i2,indx3)
				!-------------------------------------------------------------------------------
			end do
		end do
	end do
!$omp end do
!$omp end parallel



	call fftw_c2r_3D(dx_ux)
	call fftw_c2r_3D(dy_ux)
	call fftw_c2r_3D(dz_ux)

	call fftw_c2r_3D(dx_uy)
	call fftw_c2r_3D(dy_uy)
	call fftw_c2r_3D(dz_uy)

	call fftw_c2r_3D(dx_uz)
	call fftw_c2r_3D(dy_uz)
	call fftw_c2r_3D(dz_uz)

!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,dx_ux, 0)
!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,dy_ux, 0)
!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,dz_ux, 0)

!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,dx_uy, 0)
!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,dy_uy, 0)
!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,dz_uy, 0)

!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,dx_uz, 0)
!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,dy_uz, 0)
!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,dz_uz, 0)

  dx_ux = scale*dx_ux
  dy_ux = scale*dy_ux
  dz_ux = scale*dz_ux
  dx_uy = scale*dx_uy
  dy_uy = scale*dy_uy
  dz_uy = scale*dz_uy
  dx_uz = scale*dx_uz
  dy_uz = scale*dy_uz
  dz_uz = scale*dz_uz


  if (mod(cnt,20)==0) then
    ! First 'do' and 'if' are not needed if not for the 'write' statement

		do itask=0, NTask -1
			call MPI_Barrier(MPI_COMM_WORLD, ierror)
			if (itask .eq. ThisTask) then 

!!$omp parallel &
!!$omp shared(n1,n2,n3,dx_ux,dx_uy,dx_uz,dy_ux,dy_uy,dy_uz,dz_ux,dz_uy,dz_uz)&
!!$omp private(i1,i2,i3,der_mat,der_mat2,der_mat3,div,Q,R) 
!!$omp do

				do i3 = 1,local_n3,8
					do i2 = 1,n3,8
						do i1 = 1,n1,8

						der_mat(1,1) = dx_ux(i1,i2,i3)
						der_mat(1,2) = dx_uy(i1,i2,i3)
						der_mat(1,3) = dx_uz(i1,i2,i3)

						der_mat(2,1) = dy_ux(i1,i2,i3)
						der_mat(2,2) = dy_uy(i1,i2,i3)
						der_mat(2,3) = dy_uz(i1,i2,i3)

						der_mat(3,1) = dz_ux(i1,i2,i3)
						der_mat(3,2) = dz_uy(i1,i2,i3)
						der_mat(3,3) = dz_uz(i1,i2,i3)

						der_mat2 = matmul(der_mat,der_mat)
						der_mat3 = matmul(der_mat,der_mat2)  
						div = der_mat(1,1)+der_mat(2,2)+der_mat(3,3)
						Q = -0.5d0*(der_mat2(1,1)+der_mat2(2,2)+der_mat2(3,3)) 
						R = -(1.0d0/3.0d0)*(der_mat3(1,1)+der_mat3(2,2)+der_mat3(3,3))   
						!write(13,*)div,Q,R
						end do
					end do
				end do
			end if
		end do 	
!!$omp end do
!!$omp end parallel
	endif
	end subroutine dermat_serial
	
	
	
	
	
