	module mod_particle
	use mod_serial_fluid
	contains
	
	subroutine allocate_particle
	implicit none
	integer::i1

	allocate(tau(nstokes-1),one_by_tau(nstokes-1))
	open(unit = 11,file = 'stokes_no.in',status = 'old')
	print*,nstokes
	do i1 = 1, nstokes-1
		read(11,*)tau(i1)
		one_by_tau(i1) = 1.0d0/tau(i1)
	enddo
	close(11)

	allocate(Xp(nprtcl,nstokes),Yp(nprtcl,nstokes),Zp(nprtcl,nstokes))
	Xp = 0.0d0; Yp = 0.0d0; Zp = 0.0d0
	
	allocate(Uxp(nprtcl,nstokes),Uyp(nprtcl,nstokes),Uzp(nprtcl,nstokes))
	Uxp = 0.0d0; Uyp = 0.0d0; Uzp = 0.0d0 

	allocate(Vxp(nprtcl,nstokes),Vyp(nprtcl,nstokes),Vzp(nprtcl,nstokes))
	Vxp = 0.0d0; Vyp = 0.0d0; Vzp = 0.0d0 

	allocate(Axp(nprtcl,nstokes),Ayp(nprtcl,nstokes),Azp(nprtcl,nstokes))
	Axp = 0.0d0; Ayp = 0.0d0; Azp = 0.0d0 

	allocate(preAxp(nprtcl,nstokes),preAyp(nprtcl,nstokes),preAzp(nprtcl,nstokes))
	preAxp = 0.0d0; preAyp = 0.0d0; preAzp = 0.0d0 

	allocate(dAdtxp(nprtcl,nstokes),dAdtyp(nprtcl,nstokes),dAdtzp(nprtcl,nstokes))
	dAdtxp = 0.0d0; dAdtyp = 0.0d0; dAdtzp = 0.0d0 

	allocate(dx_uxp(nprtcl,nstokes),dx_uyp(nprtcl,nstokes),dx_uzp(nprtcl,nstokes))
	allocate(dy_uxp(nprtcl,nstokes),dy_uyp(nprtcl,nstokes),dy_uzp(nprtcl,nstokes))
	allocate(dz_uxp(nprtcl,nstokes),dz_uyp(nprtcl,nstokes),dz_uzp(nprtcl,nstokes))
	dx_uxp = 0.0d0; dx_uyp = 0.0d0; dx_uzp = 0.0d0 
	dy_uxp = 0.0d0; dy_uyp = 0.0d0; dy_uzp = 0.0d0 
	dz_uxp = 0.0d0; dz_uyp = 0.0d0; dz_uzp = 0.0d0 

	allocate(Qp(nprtcl,nstokes),Rp(nprtcl,nstokes),discrip(nprtcl,nstokes))
	Qp = 0.0d0; Rp = 0.0d0; discrip = 0.0d0 

	end subroutine allocate_particle
!------------------------------------------------------------------------------

	subroutine initialize_particle
	implicit none
	integer::i1,intx,inty,intz
	if (nrunpart == 1) then

	  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Vk1, 0)
		call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Vk2, 0)
		call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Vk3, 0)
!! -------normalise etc ------------------------------
	  do i1=n1+1,n1+2
		  Vk1(i1,:,:) = 0.0d0
			Vk2(i1,:,:) = 0.0d0
			Vk3(i1,:,:) = 0.0d0
	  enddo
		Vk1 = Vk1*scale
	  Vk2 = Vk2*scale
		Vk3 = Vk3*scale


	  call random_seed()
	  call random_number(Xp(:,1)); call random_number(Yp(:,1)); call random_number(Zp(:,1))

	  do i1 = 1,nprtcl
		  intx = floor(Xp(i1,1)*dfloat(n1))+1
			inty = floor(Yp(i1,1)*dfloat(n2))+1
			intz = floor(Zp(i1,1)*dfloat(n3))+1
	    Xp(i1,:) = (intx-1)*dx; Yp(i1,:) = (inty-1)*dy; Zp(i1,:) = (intz-1)*dz
		  Vxp(i1,:)=Vk1(intx,inty,intz)
			Vyp(i1,:)=Vk2(intx,inty,intz)
	    Vzp(i1,:)=Vk3(intx,inty,intz)
		enddo

	  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Vk1, 0)
		call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Vk2, 0)
		call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Vk3, 0)

 
	endif

	end subroutine initialize_particle
!------------------------------------------------------------------------

	subroutine evolv_particle 
  implicit none
  integer::ist


  call get_particle_rhs
!
!The first set of particles is always Lagrangian
!
!$OMP WORKSHARE
  Xp=Xp+Vxp*delta
  Yp=Yp+Vyp*delta
  Zp=Zp+Vzp*delta
!$OMP END WORKSHARE
  do ist=2, nstokes
!$OMP WORKSHARE
    Vxp(:,ist)=Vxp(:,ist)+Axp(:,ist)*delta
    Vyp(:,ist)=Vyp(:,ist)+Ayp(:,ist)*delta
    Vzp(:,ist)=Vzp(:,ist)+Azp(:,ist)*delta
!$OMP END WORKSHARE
  enddo

  end subroutine evolv_particle
!---------------------------------------------------------------------------

	subroutine get_particle_rhs
	implicit none
	integer::ist

	call get_velocity_at_particle_position
	
!$OMP WORKSHARE
  Vxp(:,1) = Uxp(:,1)
  Vyp(:,1) = Uyp(:,1)
  Vzp(:,1) = Uzp(:,1)
!$OMP END WORKSHARE

  do ist=2,nstokes
!$OMP WORKSHARE
    Axp(:,ist) = one_by_tau(ist-1)*(Uxp(:,ist)-Vxp(:,ist))
    Ayp(:,ist) = one_by_tau(ist-1)*(Uyp(:,ist)-Vyp(:,ist))
    Azp(:,ist) = one_by_tau(ist-1)*(Uzp(:,ist)-Vzp(:,ist))
!$OMP END WORKSHARE
  enddo

	end subroutine get_particle_rhs
!--------------------------------------------------------------------------
  subroutine get_dAdt
  implicit none
  integer::ist

  do ist=2,nstokes
!$OMP WORKSHARE
    dAdtxp(:,ist) = (Axp(:,ist)-preAxp(:,ist))/delta
    dAdtyp(:,ist) = (Ayp(:,ist)-preAyp(:,ist))/delta
    dAdtzp(:,ist) = (Azp(:,ist)-preAzp(:,ist))/delta
!$OMP END WORKSHARE
  enddo

  do ist=2,nstokes
!$OMP WORKSHARE
		preAxp(:,ist) = Axp(:,ist)
		preAyp(:,ist) = Ayp(:,ist)
		preAzp(:,ist) = Azp(:,ist)
!$OMP END WORKSHARE
  enddo

  end subroutine get_dAdt 
!--------------------------------------------------------------
	subroutine get_velocity_at_particle_position
	implicit none
	integer::ist,i1
	real*8::xxp,yyp,zzp,uuxp,uuyp,uuzp

  do ist = 1,nstokes
!$OMP PARALLEL DO PRIVATE(i1,xxp,yyp,zzp,uuxp,uuyp,uuzp) SHARED(Xp,Yp,Zp,Vk1,Vk2,Vk3,Uxp,Uyp,Uzp,n1,n2,n3,ist)
    do i1 = 1,nprtcl
      xxp = Xp(i1,ist); yyp = Yp(i1,ist); zzp = Zp(i1,ist)

      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,Vk1,uuxp)
      Uxp(i1,ist) = uuxp
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,Vk2,uuyp)
      Uyp(i1,ist) = uuyp
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,Vk3,uuzp)
      Uzp(i1,ist) = uuzp

    enddo
!$OMP END PARALLEL DO
  enddo

	end subroutine get_velocity_at_particle_position
!--------------------------------------------------------------------------

	subroutine linear_interp_to_offgrid(xpar,ypar,zpar,nn1,nn2,nn3,psi,psi_interpolated)
	implicit none
	real*8::xpar,ypar,zpar
	integer::nn1,nn2,nn3,xgrid,ygrid,zgrid,xgrid_nxt,ygrid_nxt,zgrid_nxt
	real*8,dimension(nn1+2,nn2,nn3)::psi
	real*8::psi_interpolated,a1,a2,b1,b2,c1,c2,xd,yd,zd

  call get_grid_loc(xpar,xgrid,xgrid_nxt,xd)
  call get_grid_loc(ypar,ygrid,ygrid_nxt,yd)
  call get_grid_loc(zpar,zgrid,zgrid_nxt,zd)
	
	a1=(1.0d0-zd)*psi(xgrid,ygrid,zgrid)+zd*psi(xgrid,ygrid,zgrid_nxt)
	a2=(1.0d0-zd)*psi(xgrid,ygrid_nxt,zgrid)+zd*psi(xgrid,ygrid_nxt,zgrid_nxt)
	b1=(1.0d0-zd)*psi(xgrid_nxt,ygrid,zgrid)+zd*psi(xgrid_nxt,ygrid,zgrid_nxt)
  b2=(1.0d0-zd)*psi(xgrid_nxt,ygrid_nxt,zgrid)+zd*psi(xgrid_nxt,ygrid_nxt,zgrid_nxt)

  c1=(1.0d0-yd)*a1+yd*a2
  c2=(1.0d0-yd)*b1+yd*b2

  psi_interpolated=(1.0d0-xd)*c1+xd*c2

	end subroutine linear_interp_to_offgrid
!--------------------------------------------------------------------------------------------------------------------------

	subroutine get_grid_loc(xpar,grid,grid_nxt,xd)
	implicit none
	real*8::xpar,xd
	integer::grid,grid_nxt
	real*8::Xp_fld

	Xp_fld=modulo(xpar,(2.0d0*pi))	
	grid=floor(Xp_fld/dx)+1
	grid_nxt=(1-grid/n1)*grid+1
	xd=modulo(Xp_fld,dx)/dx

	end subroutine get_grid_loc
!---------------------------------------------------------------------------------------------------------------------

	subroutine dermat_at_particle
	implicit none
	integer::ist,i1
	real*8::interpolated_value,xxp,yyp,zzp
	real*8,dimension(3,3)::DM,DM_sqr,DM_cub
	
  do ist = 1,nstokes
!$OMP PARALLEL DO PRIVATE(i1,xxp,yyp,zzp,interpolated_value,DM,DM_sqr,DM_cub) &
!$OMP SHARED(Xp,Yp,Zp,dx_ux,dx_uy,dx_uz,dy_ux,dy_uy,dy_uz,dz_ux,dz_uy,dz_uz,dx_uxp,dx_uyp, &
!$OMP dx_uzp,dy_uxp,dy_uyp,dy_uzp,dz_uxp,dz_uyp,dz_uzp,n1,n2,n3,ist,Qp,Rp,discrip)

    do i1 = 1,nprtcl
      xxp = Xp(i1,ist); yyp = Yp(i1,ist); zzp = Zp(i1,ist)

      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dx_ux,interpolated_value)
      dx_uxp(i1,ist) = interpolated_value
			DM(1,1) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dx_uy,interpolated_value)
      dx_uyp(i1,ist) = interpolated_value
			DM(1,2) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dx_uz,interpolated_value)
      dx_uzp(i1,ist) = interpolated_value
			DM(1,3) = interpolated_value

      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dy_ux,interpolated_value)
      dy_uxp(i1,ist) = interpolated_value
			DM(2,1) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dy_uy,interpolated_value)
      dy_uyp(i1,ist) = interpolated_value
			DM(2,2) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dy_uz,interpolated_value)
      dy_uzp(i1,ist) = interpolated_value
			DM(2,3) = interpolated_value

      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dz_ux,interpolated_value)
      dz_uxp(i1,ist) = interpolated_value
			DM(3,1) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dz_uy,interpolated_value)
      dz_uyp(i1,ist) = interpolated_value
			DM(3,2) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dz_uz,interpolated_value)
      dz_uzp(i1,ist) = interpolated_value
			DM(3,3) = interpolated_value
			
			DM_sqr = matmul(DM,DM)
			DM_cub = matmul(DM,DM_sqr)
			
			Qp(i1,ist) = -0.5d0*(DM_sqr(1,1)+DM_sqr(2,2)+DM_sqr(3,3))
			Rp(i1,ist) = -(1.0d0/3.0d0)*(DM_cub(1,1)+DM_cub(2,2)+DM_cub(3,3))
			discrip(i1,ist) = Qp(i1,ist)**3+(27.0d0/4.0d0)*Rp(i1,ist)**2

    enddo
!$OMP END PARALLEL DO
  enddo
	
	end subroutine dermat_at_particle
!-----------------------------------------------------------------------------------------------------------

	subroutine open_output_particle
	implicit none
  character*80 :: fname1,fname2
  integer::i1,fno,ist
!!
  do ist = 1,nstokes
    write(fname1,'(g8.0)')  ist
    fno = 1000*ist + 30
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/xptrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/yptrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/zptrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/vxint.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/vyint.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/vzint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/uxint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/uyint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/uzint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dAdtxint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dAdtyint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dAdtzint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/R_int.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/Delta_int.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dx_uxp.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dx_uyp.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dx_uzp.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dy_uxp.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dy_uyp.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dy_uzp.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dz_uxp.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dz_uyp.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dz_uzp.out',status='unknown')
  enddo
!!
	end subroutine open_output_particle
!-------------------------------------------------------------------------------------------------

	subroutine write_output_particle
	implicit none
	integer::i1,fno,ist
	character*500 :: cnpar,formp

  write(cnpar,'(g8.0)') nprtcl
  formp = '('//trim(adjustl(cnpar))//'es16.3E2)' 
	
  do ist = 1,nstokes
    fno = 1000*ist + 30
    write(fno,formp) (Xp(i1,ist),i1=1,nprtcl);  fno = fno + 1
    write(fno,formp) (Yp(i1,ist),i1=1,nprtcl);  fno = fno + 1
    write(fno,formp) (Zp(i1,ist),i1=1,nprtcl);  fno = fno + 1
    write(fno,formp) (Vxp(i1,ist),i1=1,nprtcl); fno = fno + 1
    write(fno,formp) (Vyp(i1,ist),i1=1,nprtcl); fno = fno + 1
    write(fno,formp) (Vzp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (Uxp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (Uyp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (Uzp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dAdtxp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dAdtyp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dAdtzp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (Rp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (discrip(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dx_uxp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dx_uyp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dx_uzp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dy_uxp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dy_uyp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dy_uzp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dz_uxp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dz_uyp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dz_uzp(i1,ist),i1=1,nprtcl) 
  enddo

		
	end subroutine write_output_particle

	end module mod_particle
