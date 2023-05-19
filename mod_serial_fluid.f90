!!             mod_serial_fluid.f90
!! This is the module containing the defination and commons for the 
!! serial code. As all the serial ffts requal arrays of same type 
!! we expect this should not be changed as we go from one machine  
!! to another.
	module mod_serial_fluid
	use, intrinsic :: iso_c_binding		
	implicit none
	save
!! ----velocites  -------------
	real*8,allocatable,dimension(:,:,:) :: Vk1,Vk2,Vk3,VWk1,VWk2,& 
                                        VWk3,Wk1,Wk2,Wk3,tm1,tm2,tm3,&
                                        fk1,fk2,fk3
                                                 
                                                                  
!----------lagrangian trajectaries-----------------------------------
	real*8,allocatable,dimension(:,:)::Xp,Yp,Zp,Uxp,Uyp,Uzp,Vxp,Vyp,Vzp,Axp,Ayp,Azp,Qp,Rp,discrip, &
																			dx_uxp,dx_uyp,dx_uzp,dy_uxp,dy_uyp,dy_uzp,dz_uxp,dz_uyp,dz_uzp, &
																			preAxp,preAyp,preAzp,dAdtxp,dAdtyp,dAdtzp
	real*8,allocatable,dimension(:)::tau,one_by_tau
  logical::forcing,particl
  real*8,allocatable,dimension(:,:,:)::dx_ux,dx_uy,dx_uz, &
                              dy_ux,dy_uy,dy_uz,dz_ux,dz_uy,dz_uz

	real*8,allocatable,dimension(:)::uxe,uye,uze,dx_uxe,dx_uye,dx_uze,dy_uxe,dy_uye,dy_uze,dz_uxe,dz_uye,dz_uze

	integer,allocatable,dimension(:,:,:):: mask	
	integer,allocatable,dimension(:):: N_modes	
	real*8:: D,amp
!! --energy,viscosity,forced modes etc --------------
	real*8,allocatable,dimension(:,:) :: E_Omega, E_Omega_total  
	double precision,allocatable,dimension(:) :: iniEk 
	integer,allocatable,dimension(:,:) :: forced_modes1,forced_modes2,forced_modes3
	integer :: count1,count2,count3
!! ---parameters of simulations -----------------
	real*8 :: edr,tlr_micro_scale,dissipation_scale,vkrms,Rlambda, &
              tau_eddy,integral_scale,dx,dy,dz,t1,t2,t, t_all, dt_by_dx, kini, fe3
!! -----real space structure functions -----------------------
	real*8,allocatable,dimension(:,:) :: S_p
	integer :: pmax,nthreads,count_out,cnt
!! -----common parameters ----------------------------
	integer :: thousand,nn,maxiter,nrun,n1d, &
               nshell,navg,nalias_sqr,nprtcl,nstokes,nrunpart,vel_writer
	real*8 :: constA,constB,vis,vis2,delta,pi,length,factor,umax1,umax2,umax3,umax,c_cfl
	double precision :: fixed_energy1,fixed_energy2,fixed_energy3
	complex*16 :: zi,zone 
!! -----global arrays --------------
	real*8,allocatable,dimension(:) :: time_increment, time_increment_total
	integer,allocatable,dimension(:) :: den_state, den_state_total
!! ----------------------------------------------------------


!process related
integer :: ThisTask, NTask, ierror, itask, iret, n_proc
character(100) smallstring, cur_dir
character(200) bigstring 

!data related
integer (C_INTPTR_T) :: n1,n2,n3, n1h,n1hf
integer (C_INTPTR_T) :: alloc_local, local_n3_offset, local_n3
real(C_DOUBLE), pointer :: r_data(:,:,:)
complex(C_DOUBLE_COMPLEX), pointer :: c_data(:,:,:)
integer*8 :: n3_low, n3_high, indx3

integer :: data_scatter

!! -----------For FFTW --------------------------------------
	type(C_PTR) :: plan_r2c, plan_c2r, cdatar, cdatac
	integer,dimension(3) :: dim
	double precision :: scale
	include 'mpif.h'
	include "fftw3-mpi.f03"


	
!! --------------------------------------------------------
			end module mod_serial_fluid 
