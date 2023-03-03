!!      serial_spectra.f90
!! Calculates the spectra for a given velocity field, this is the
!! serial version. This should not change as we move from one 
!! serial machine to another
	subroutine energy_serial 
	use mod_serial_fluid
	implicit none
	real*8 :: energy_norm,rk2,rk,energy,dissipation,one_by_ncube, &
				volume_element	
	integer :: i1,i2,i3,k1,k2,k3,ksqr,mshl,ireal,iimag 
!! --------------------------------------------------------------
	one_by_ncube = 1.0d0/dfloat(n1*n2*n3)
	volume_element = one_by_ncube*(length**3.0d0)
	energy_norm = one_by_ncube*one_by_ncube 

!$omp parallel &
!$omp shared(Vk1,Vk2,Vk3,n2,n3,n1hf,energy_norm) &
!$omp private(i1,i2,i3,k1,k2,k3,indx3,ksqr,rk2,rk,mshl,ireal,iimag,energy,dissipation) &
!$omp reduction(+ : E_Omega)
!$omp do

	do i3 = n3_low,n3_high
		indx3 = i3 - n3_low + 1
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
				energy = energy_norm*( &
			Vk1(ireal,i2,indx3)**2 + Vk1(iimag,i2,indx3)**2 + &
			Vk2(ireal,i2,indx3)**2 + Vk2(iimag,i2,indx3)**2 + &
			Vk3(ireal,i2,indx3)**2 + Vk3(iimag,i2,indx3)**2 )
				dissipation = energy*rk2
			if((k1.eq.0).or.(k1.eq.n1h)) then
				E_Omega(mshl,1) = E_Omega(mshl,1) + energy
				E_Omega(mshl,2) = E_Omega(mshl,2) + dissipation
			else
				E_Omega(mshl,1) = E_Omega(mshl,1) + 2.0d0*energy
				E_Omega(mshl,2) = E_Omega(mshl,2) + 2.0d0*dissipation
				endif	
!! ----------------- 
			enddo
		enddo
	enddo

!$omp end do
!$omp end parallel
			call MPI_Allreduce(E_Omega , E_Omega_total, (nshell+1)*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierror)
			E_Omega = 0.5d0*E_Omega_total

!! --------------------------------------------------------
			end subroutine energy_serial 
