!*******************************************************************
!*   *     Subroutine to calculate energy injection rate    *     *!
!*******************************************************************

  subroutine energy_injection
	use mod_serial_fluid
	use mod_initial
	implicit none

	integer::i1,i2,i3
	real*8::f_dot_v,sum_f_dot_v, sum_f_dot_v_total
	
	call fftw_c2r_3D(fk1)
	call fftw_c2r_3D(fk2)
	call fftw_c2r_3D(fk3)

!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,fk1, 0)
!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,fk2, 0)
!  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,fk3, 0)

	 do i1=n1+1,n1+2
		fk1(i1,:,:) = 0.0d0
		fk2(i1,:,:) = 0.0d0
		fk3(i1,:,:) = 0.0d0
	enddo
	fk1 = fk1*scale
	fk2 = fk2*scale
	fk3 = fk3*scale

	sum_f_dot_v = 0.0d0
!$omp parallel &
!$omp shared(fk1,fk2,fk3,Vk1,Vk2,Vk3,n1,n2,n3) &
!$omp private(i1,i2,i3,f_dot_v) &
!$omp reduction(+ : sum_f_dot_v)
!$omp do

	do i3 = 1,local_n3
		do i2 = 1,n2
			do i1 = 1,n1

			f_dot_v = fk1(i1,i2,i3)*Vk1(i1,i2,i3)+ &
                                  fk2(i1,i2,i3)*Vk2(i1,i2,i3)+ &
                                  fk3(i1,i2,i3)*Vk3(i1,i2,i3)

			sum_f_dot_v = sum_f_dot_v+f_dot_v

			enddo
		enddo
	enddo
!$omp end do
!$omp end parallel

	call MPI_Allreduce(sum_f_dot_v , sum_f_dot_v_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierror)
	sum_f_dot_v = sum_f_dot_v_total*scale
	if (ThisTask .eq. 0) then 
		write(202,*)sum_f_dot_v
	end if 
  end subroutine energy_injection
  
  
  
  
  
  
  
  
  
  
  
