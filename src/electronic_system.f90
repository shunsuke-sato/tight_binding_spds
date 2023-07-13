module electronic_system
  use mpi
  use parallel
  use communication
  use math
  use constants
  use inputoutput
  implicit none
  private

  public :: initialize_electronic_system, &
            calc_bandstructure_zincblende, &
            set_equilibrium_density_matrix, &
            dt_evolve_elec_system, &
            dt_evolve_elec_system_mod, &
            calc_current, &
            calc_num_electron

! material class
  integer,parameter :: n_diamond = 0, n_zincblende = 1
  integer,parameter :: material_class = n_zincblende

! tight-binding
  real(8) :: lattice_vec(3,3), lattice_const, volume
  real(8) :: reciprocal_lattice_vec(3,3)
  integer,parameter :: nspin = 2
  integer :: nband, nkpoint
  integer :: nk1, nk2, nk3
  integer :: nk_s, nk_e
  real(8),allocatable :: kvec0(:,:),kvec(:,:)
  complex(8),allocatable :: zpsi(:,:),zrho_dm(:,:,:)
  complex(8),allocatable :: zHam_mat(:,:,:)
  complex(8),allocatable :: zJ_mat(:,:,:,:)

  integer :: num_nearest_neighbor
  real(8),allocatable :: Rvec_ac(:,:) !! vector from anion to cation


!tight binding parameters zinc blende
  real(8) :: energy_unit_tb
  real(8) :: Ea_s,Ea_p,Ea_d,Ea_s2
  real(8) :: Ec_s,Ec_p,Ec_d,Ec_s2

  real(8) :: ss_sigma
  real(8) :: s2s2_sigma
  real(8) :: s2a_sc_sigma
  real(8) :: sa_s2c_sigma

  real(8) :: sapc_sigma
  real(8) :: scpa_sigma
  real(8) :: s2apc_sigma
  real(8) :: s2cpa_sigma

  real(8) :: sadc_sigma
  real(8) :: scda_sigma
  real(8) :: s2adc_sigma
  real(8) :: s2cda_sigma

  real(8) :: pp_sigma
  real(8) :: pp_pi

  real(8) :: padc_sigma
  real(8) :: pcda_sigma
  real(8) :: padc_pi
  real(8) :: pcda_pi

  real(8) :: dd_sigma
  real(8) :: dd_pi
  real(8) :: dd_delta

  real(8) :: delta_a
  real(8) :: delta_c

! two-body integral
  real(8),allocatable :: E2c_int(:,:,:)


! relaxation
  real(8) :: T1_relax, T2_relax

! focal spot average
  logical :: if_focal_spot_average 
  integer :: num_sample_focal_spot_average
  real(8),allocatable :: alpha_fsa(:)
  integer :: npower_focal_spot_average

! band freezing
  logical :: if_band_frozen
  real(8),allocatable :: blocking_matrix_band_frozen(:,:)

contains
!----------------------------------------------------------------------------
subroutine initialize_electronic_system
  implicit none
  real(8) :: vec_t(3)
  integer :: nk_average, nk_remainder
  integer :: ik1,ik2,ik3,ik
  real(8) :: T1_relax_fs, T2_relax_fs

  nband = 2*2*(1+3+5+1) ! 40 bands
  lattice_const = 5.4635d0*angstrom


  call read_basic_input('nk1',nk1,val_default = -1)
  call read_basic_input('nk2',nk2,val_default = -1)
  call read_basic_input('nk3',nk3,val_default = -1)

  nkpoint = nk1*nk2*nk3

  call read_basic_input('if_focal_spot_average',if_focal_spot_average,val_default = .false.)
  call read_basic_input('num_sample_focal_spot_average',num_sample_focal_spot_average,val_default = -1)
  call read_basic_input('npower_focal_spot_average',npower_focal_spot_average,val_default = 1)


  if(if_focal_spot_average)then
    nkpoint = 12*num_sample_focal_spot_average
  end if

  call read_basic_input('if_band_frozen',if_band_frozen,val_default = .false.)
  allocate(blocking_matrix_band_frozen(nband,nband))
  blocking_matrix_band_frozen(:,:) = 1d0
  if(if_band_frozen)call set_blocking_matrix_band_frozen

  nk_average = nkpoint/comm_nproc_global
  nk_remainder = mod(nkpoint,comm_nproc_global)
  if(comm_id_global+1 <= nk_remainder)then
    nk_s = 1 + comm_id_global*(nk_average+1)
    nk_e   = nk_s + (nk_average + 1) -1
  else
    nk_s = 1 + nk_remainder*(nk_average+1) + nk_average*(comm_id_global - nk_remainder)
    nk_e    = nk_s + nk_average  -1
  end if


  allocate(kvec0(3,nk_s:nk_e),kvec(3,nk_s:nk_e),alpha_fsa(nk_s:nk_e))
  allocate(zpsi(nband,nk_s:nk_e), zHam_mat(nband,nband,nk_s:nk_e))
  allocate(zJ_mat(nband,nband,3,nk_s:nk_e))

! nband = 1-20
!! anion spin-up
! 1 = s
! 2 = px
! 3 = py
! 4 = pz
! 5 = dxy
! 6 = dyz
! 7 = dzx
! 8 = dx2-y2
! 9 = d3z2-r2
!10 = s*
!! cation spin-up
!11 = s
!12 = px
!13 = py
!14 = pz
!15 = dxy
!16 = dyz
!17 = dzx
!18 = dx2-y2
!19 = d3z2-r2
!20 = s*
!! spin down components (21-40)

  lattice_vec(1,1) = 0d0
  lattice_vec(2,1) = 1d0
  lattice_vec(3,1) = 1d0

  lattice_vec(1,2) = 1d0
  lattice_vec(2,2) = 0d0
  lattice_vec(3,2) = 1d0

  lattice_vec(1,3) = 1d0
  lattice_vec(2,3) = 1d0
  lattice_vec(3,3) = 0d0
  
  lattice_vec = 0.5d0*lattice_const*lattice_vec
  vec_t = cross_product(lattice_vec(:,2),lattice_vec(:,3))
  volume = abs(sum(lattice_vec(:,1)*vec_t(:)))
!  write(*,*)"volume=",volume

  reciprocal_lattice_vec(:,1)=2d0*pi/volume*cross_product(lattice_vec(:,2),lattice_vec(:,3))
  reciprocal_lattice_vec(:,2)=2d0*pi/volume*cross_product(lattice_vec(:,3),lattice_vec(:,1))
  reciprocal_lattice_vec(:,3)=2d0*pi/volume*cross_product(lattice_vec(:,1),lattice_vec(:,2))

  num_nearest_neighbor = 4
  allocate(Rvec_ac(3,num_nearest_neighbor))
  allocate(E2c_int(10,10,num_nearest_neighbor))

!  Rvec_ac(:,1) = lattice_const*0.25d0
  Rvec_ac(:,1) = lattice_const*0.25d0
  Rvec_ac(:,2) = Rvec_ac(:,1) - lattice_vec(:,1)
  Rvec_ac(:,3) = Rvec_ac(:,1) - lattice_vec(:,2)
  Rvec_ac(:,4) = Rvec_ac(:,1) - lattice_vec(:,3)

!  include "include_tb_parameters/set_GaAs_Jancu1998.f90"
  include "include_tb_parameters/set_GaAs_Tan2013.f90"

  if(if_focal_spot_average)then
    call prepare_sampling_for_focal_spot_average
  else ! uniform sampling
    alpha_fsa = 1d0
    kvec0 = 0d0; kvec=0d0
    ik = 0
    do ik1 = 1, nk1
      do ik2 = 1, nk2
        do ik3 = 1, nk3
          ik = ik + 1
          if(ik >= nk_s .and. ik <= nk_e)then
            kvec0(:,ik) = (ik1-1)/dble(nk1)*reciprocal_lattice_vec(:,1) &
                +(ik2-1)/dble(nk2)*reciprocal_lattice_vec(:,2) &
                +(ik3-1)/dble(nk3)*reciprocal_lattice_vec(:,3) 
          end if
        end do
      end do
    end do
    kvec = kvec0
  end if


  call read_basic_input('T1_relax_fs',T1_relax_fs,val_default = -1d0)
  call read_basic_input('T2_relax_fs',T2_relax_fs,val_default = -1d0)
  T1_relax = T1_relax_fs*fs
  T2_relax = T2_relax_fs*fs

contains
  function cross_product(vec1, vec2)
    real(8) :: cross_product(3)
    real(8) :: vec1(3), vec2(3)

    cross_product(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
    cross_product(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
    cross_product(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)

  end function cross_product

end subroutine initialize_electronic_system
!----------------------------------------------------------------------------
subroutine calc_two_center_integral
  implicit none
  real(8) :: l,m,n
  integer :: i

  E2c_int = 0d0

  do i = 1, num_nearest_neighbor
    l = Rvec_ac(1,i)/sqrt(sum(Rvec_ac(:,i)**2))
    m = Rvec_ac(2,i)/sqrt(sum(Rvec_ac(:,i)**2))
    n = Rvec_ac(3,i)/sqrt(sum(Rvec_ac(:,i)**2))

!    write(*,*)l,m,n

! anion s-orbital
    E2c_int(1,1,i) = ss_sigma !ss
    E2c_int(1,2,i) = l*sapc_sigma !sx
    E2c_int(1,3,i) = m*sapc_sigma !sy
    E2c_int(1,4,i) = n*sapc_sigma !sz

    E2c_int(1,5,i) = sqrt(3d0)*l*m*sadc_sigma !s,xy
    E2c_int(1,6,i) = sqrt(3d0)*m*n*sadc_sigma !s,yx
    E2c_int(1,7,i) = sqrt(3d0)*n*l*sadc_sigma !s,zx
    E2c_int(1,8,i) = 0.5d0*sqrt(3d0)*(l**2-m**2)*sadc_sigma !s,x2-y2
    E2c_int(1,9,i) = (n**2-0.5d0*(l**2+m**2))*sadc_sigma !s,3z2-r2

    E2c_int(1,10,i) = sa_s2c_sigma !ss2

! anion px-orbital
    E2c_int(2,1,i) = (-1d0)*l*scpa_sigma !x,s

    E2c_int(2,2,i) = l**2*pp_sigma+(1d0-l**2)*pp_pi !x,x
    E2c_int(2,3,i) = l*m*pp_sigma-l*m*pp_pi !x,y
    E2c_int(2,4,i) = l*n*pp_sigma-l*n*pp_pi !x,z

    E2c_int(2,5,i)  = sqrt(3d0)*l**2*m*padc_sigma+m*(1d0-2d0*l**2)*padc_pi !x,xy
    E2c_int(2,6,i)  = sqrt(3d0)*l*m*n*padc_sigma-2d0*l*m*n*padc_pi !x, yz
    E2c_int(2,7,i)  = sqrt(3d0)*l**2*n*padc_sigma+n*(1d0-2d0*l**2)*padc_pi !x, zx

    E2c_int(2,8,i)  = 0.5d0*sqrt(3d0)*l*(l**2-m**2)*padc_sigma &
                    +l*(1d0-l**2+m**2)*padc_pi !x,x2-y2
    E2c_int(2,9,i)  = l*(n**2-0.5d0*(l**2+m**2))*padc_sigma &
        -sqrt(3d0)*l*n**2*padc_pi !x, 3z2-r2
    E2c_int(2,10,i) = (-1d0)*l*s2cpa_sigma !x, s2

! anion py-orbital
    E2c_int(3,1,i) = (-1d0)*(m*scpa_sigma) !y,s
    E2c_int(3,2,i) = l*m*pp_sigma-l*m*pp_pi !y,x
    E2c_int(3,3,i) = m**2*pp_sigma+(1d0-m**2)*pp_pi !y,y
    E2c_int(3,4,i) = m*n*pp_sigma-m*n*pp_pi !y,z
    E2c_int(3,5,i) = sqrt(3d0)*m**2*l*padc_sigma+l*(1d0-2d0*m**2)*padc_pi !y,xy
    E2c_int(3,6,i) = sqrt(3d0)*m**2*n*padc_sigma+n*(1d0-2d0*m**2)*padc_pi !y,yz
    E2c_int(3,7,i) = sqrt(3d0)*m*n*l*padc_sigma-2d0*m*n*l*padc_pi !y,zx
    E2c_int(3,8,i) = 0.5d0*sqrt(3d0)*m*(l**2-m**2)*padc_sigma &
                     -m*(1d0+l**2-m**2)*padc_pi !y,x2-y2
    E2c_int(3,9,i) = m*(n**2-0.5d0*(l**2+m**2))*padc_sigma &
                     -sqrt(3d0)*m*n**2*padc_pi !y,3z2-r2
    E2c_int(3,10,i) = (-1d0)*m*s2cpa_sigma !y,s2

! anion pz-orbital
    E2c_int(4,1,i)  = (-1d0)*n*scpa_sigma !z,s
    E2c_int(4,2,i)  = n*l*pp_sigma-n*l*pp_pi !z,x
    E2c_int(4,3,i)  = m*n*pp_sigma-m*n*pp_pi !z,y
    E2c_int(4,4,i)  = n**2*pp_sigma+(1d0-n**2)*pp_pi !z,z
    E2c_int(4,5,i)  = sqrt(3d0)*l*m*n*padc_sigma-2d0*l*m*n*padc_pi !z,xy
    E2c_int(4,6,i)  = sqrt(3d0)*n**2*m*padc_sigma+m*(1d0-2d0*n**2)*padc_pi !z,yz
    E2c_int(4,7,i)  = sqrt(3d0)*n**2*l*padc_sigma+l*(1d0-2d0*n**2)*padc_pi !z,zx
    E2c_int(4,8,i)  = 0.5d0*sqrt(3d0)*n*(l**2-m**2)*padc_sigma &
                      -n*(l**2-m**2)*padc_pi !z,x2-y2
    E2c_int(4,9,i)  = n*(n**2-0.5d0*(l**2+m**2))*padc_sigma &
                      +sqrt(3d0)*n*(l**2+m**2)*padc_pi !z,3z2-r2
    E2c_int(4,10,i) = (-1d0)*n*s2cpa_sigma !z,s2

! anion dxy-orbital
    E2c_int(5,1,i)  =  sqrt(3d0)*l*m*scda_sigma !xy,s
    E2c_int(5,2,i)  =  (-1d0)*(sqrt(3d0)*l**2*m*pcda_sigma &
                       +m*(1d0-2d0*l**2)*pcda_pi) !xy,x
    E2c_int(5,3,i)  =  (-1d0)*(sqrt(3d0)*m**2*l*pcda_sigma &
                       +l*(1d0-2d0*m**2)*pcda_pi) !xy,y
    E2c_int(5,4,i)  =  (-1d0)*(sqrt(3d0)*l*m*n*pcda_sigma &
                       -2d0*l*m*n*pcda_pi) !xy,z
    E2c_int(5,5,i)  =  3d0*l**2*m**2*dd_sigma+(l**2+m**2-4d0*l**2*m**2)*dd_pi &
                       +(n**2+l**2*m**2)*dd_delta !xy,xy
    E2c_int(5,6,i)  =  3d0*l*m**2*n*dd_sigma+l*n*(1d0-4d0*m**2)*dd_pi &
                       +l*n*(m**2-1d0)*dd_delta !xy,yz
    E2c_int(5,7,i)  =  3d0*l**2*m*n*dd_sigma+m*n*(1d0-4d0*l**2)*dd_pi &
                       +m*n*(l**2-1d0)*dd_delta !xy,zx
    E2c_int(5,8,i)  =  1.5d0*l*m*(l**2-m**2)*dd_sigma+2d0*l*m*(m**2-l**2)*dd_pi &
                       +0.5d0*l*m*(l**2-m**2)*dd_delta !xy,x2-y2
    E2c_int(5,9,i)  =  sqrt(3d0)*l*m*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                       -2d0*sqrt(3d0)*l*m*n**2*dd_pi &
                       +0.5d0*sqrt(3d0)*l*m*(1d0+n**2)*dd_delta !xy,3z2-r2
    E2c_int(5,10,i) =  sqrt(3d0)*l*m*s2cda_sigma !xy,s2

! anion dyz-orbital
    E2c_int(6,1,i)  = sqrt(3d0)*m*n*scda_sigma !yz,s
    E2c_int(6,2,i)  = (-1d0)*(sqrt(3d0)*l*m*n*pcda_sigma &
                      -2d0*l*m*n*pcda_pi) !yz,x
    E2c_int(6,3,i)  = (-1d0)*(sqrt(3d0)*m**2*n*pcda_sigma &
                      +n*(1d0-2d0*m**2)*pcda_pi) !yz,y
    E2c_int(6,4,i)  = (-1d0)*(sqrt(3d0)*n**2*m*pcda_sigma &
                      +m*(1d0-2d0*n**2)*pcda_pi) !yz,z
    E2c_int(6,5,i)  = 3d0*l*m**2*n*dd_sigma+l*n*(1d0-4d0*m**2)*dd_pi &
                      +l*n*(m**2-1d0)*dd_delta !yz,xy
    E2c_int(6,6,i)  = 3d0*m**2*n**2*dd_sigma+(m**2+n**2-4d0*m**2*n**2)*dd_pi &
                      +(l**2+m**2*n**2)*dd_delta !yz,yz
    E2c_int(6,7,i)  = 3d0*m*n**2*l*dd_sigma+m*l*(1d0-4d0*n**2)*dd_pi &
                      +m*l*(n**2-1d0)*dd_delta !yz,zx
    E2c_int(6,8,i)  = 1.5d0*m*n*(l**2-m**2)*dd_sigma &
                      -m*n*(1d0+2d0*(l**2-m**2))*dd_pi &
                      +m*n*(1d0+0.5d0*(l**2-m**2))*dd_delta !yz,x2-y2
    E2c_int(6,9,i)  = sqrt(3d0)*m*n*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                      +sqrt(3d0)*m*n*(l**2+m**2-n**2)*dd_pi &
                      -0.5d0*sqrt(3d0)*m*n*(l**2+m**2)*dd_delta !yz,3z2-r2
    E2c_int(6,10,i) = sqrt(3d0)*m*n*s2cda_sigma !yz,s2

! anion dzx-orbital
    E2c_int(7,1,i)  = sqrt(3d0)*n*l*scda_sigma  !zx,s
    E2c_int(7,2,i)  = (-1d0)*(sqrt(3d0)*l**2*n*pcda_sigma &
                      +n*(1d0-2d0*l**2)*pcda_pi) !zx,x
    E2c_int(7,3,i)  = (-1d0)*(sqrt(3d0)*l*m*n*pcda_sigma &
                      -2d0*l*m*n*pcda_pi) !zx,y
    E2c_int(7,4,i)  = (-1d0)*(sqrt(3d0)*n**2*l*pcda_sigma &
                      +l*(1d0-2d0*n**2)*pcda_pi) !zx,z
    E2c_int(7,5,i)  = 3d0*l**2*m*n*dd_sigma+m*n*(1d0-4d0*l**2)*dd_pi &
                      +m*n*(l**2-1d0)*dd_delta !zx,xy
    E2c_int(7,6,i)  = 3d0*n**2*l*m*dd_sigma+l*m*(1d0-4d0*n**2)*dd_pi &
                      +l*m*(n**2-1d0)*dd_delta !zx,yz
    E2c_int(7,7,i)  = 3d0*n**2*l**2*dd_sigma+(n**2+l**2-4d0*n**2*l**2)*dd_pi &
                      +(m**2+n**2*l**2)*dd_delta !zx,zx
    E2c_int(7,8,i)  = 1.5d0*n*l*(l**2-m**2)*dd_sigma &
                      +n*l*(1d0-2d0*(l**2-m**2))*dd_pi &
                      -n*l*(1d0-0.5d0*(l**2-m**2))*dd_delta !zx,x2-y2
    E2c_int(7,9,i)  = sqrt(3d0)*l*n*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                      +sqrt(3d0)*l*n*(l**2+m**2-n**2)*dd_pi &
                      -0.5d0*sqrt(3d0)*l*n*(l**2+m**2)*dd_delta !zx,3z2-r2
    E2c_int(7,10,i) = sqrt(3d0)*n*l*s2cda_sigma !zx,s2

! anion dx2-y2-orbital
    E2c_int(8,1,i)  = 0.5d0*sqrt(3d0)*(l**2-m**2)*scda_sigma !x2-y2,s
    E2c_int(8,2,i)  = (-1d0)*(0.5d0*sqrt(3d0)*l*(l**2-m**2)*pcda_sigma &
                      +l*(1d0-l**2+m**2)*pcda_pi) !x2-y2,x
    E2c_int(8,3,i)  = (-1d0)*(0.5d0*sqrt(3d0)*m*(l**2-m**2)*pcda_sigma &
                      -m*(1d0+l**2-m**2)*pcda_pi) !x2-y2,y
    E2c_int(8,4,i)  = (-1d0)*(0.5d0*sqrt(3d0)*n*(l**2-m**2)*pcda_sigma &
                      -n*(l**2-m**2)*pcda_pi) !x2-y2,z
    E2c_int(8,5,i)  = 1.5d0*l*m*(l**2-m**2)*dd_sigma+2d0*l*m*(m**2-l**2)*dd_pi &
                      +0.5d0*l*m*(l**2-m**2)*dd_delta !x2-y2,xy
    E2c_int(8,6,i)  = 1.5d0*m*n*(l**2-m**2)*dd_sigma-m*n*(1d0+2d0*(l**2-m**2))*dd_pi &
                      +m*n*(1d0+0.5d0*(l**2-m**2))*dd_delta !x2-y2,yz
    E2c_int(8,7,i)  = 1.5d0*n*l*(l**2-m**2)*dd_sigma+n*l*(1d0-2d0*(l**2-m**2))*dd_pi &
                      -n*l*(1d0-0.5d0*(l**2-m**2))*dd_delta !x2-y2,zx
    E2c_int(8,8,i)  = 0.75d0*(l**2-m**2)**2*dd_sigma &
                      +(l**2+m**2-(l**2-m**2)**2)*dd_pi &
                      +(n**2+0.25d0*(l**2-m**2)**2)*dd_delta !x2-y2,x2-y2
    E2c_int(8,9,i)  = 0.5d0*sqrt(3d0)*(l**2-m**2)*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                      +sqrt(3d0)*n**2*(m**2-l**2)*dd_pi &
                      +0.25d0*sqrt(3d0)*(1d0+n**2)*(l**2-m**2)*dd_delta !x2-y2,3z2-r2
    E2c_int(8,10,i) = 0.5d0*sqrt(3d0)*(l**2-m**2)*s2cda_sigma !x2-y2,s2

! anion d3z2-r2-orbital
    E2c_int(9,1,i)  = (n**2-0.5d0*(l**2+m**2))*scda_sigma !3z2-r2,s
    E2c_int(9,2,i)  = (-1d0)*(l*(n**2-0.5d0*(l**2+m**2))*pcda_sigma &
                      -sqrt(3d0)*l*n**2*pcda_pi) !3z2-r2,x
    E2c_int(9,3,i)  = (-1d0)*(m*(n**2-0.5d0*(l**2+m**2))*pcda_sigma &
                      -sqrt(3d0)*m*n**2*pcda_pi) !3z2-r2,y
    E2c_int(9,4,i)  = (-1d0)*(n*(n**2-0.5d0*(l**2+m**2))*pcda_sigma &
                      +sqrt(3d0)*n*(l**2+m**2)*pcda_pi) !3z2-r2,z
    E2c_int(9,5,i)  = sqrt(3d0)*l*m*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                      -2d0*sqrt(3d0)*l*m*n**2*dd_pi &
                      +0.5d0*sqrt(3d0)*l*m*(1d0+n**2)*dd_delta !3z2-r2,xy
    E2c_int(9,6,i)  = sqrt(3d0)*m*n*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                      +sqrt(3d0)*m*n*(l**2+m**2-n**2)*dd_pi &
                      -0.5d0*sqrt(3d0)*m*n*(l**2+m**2)*dd_delta !3z2-r2,yz
    E2c_int(9,7,i)  = sqrt(3d0)*l*n*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                      +sqrt(3d0)*l*n*(l**2+m**2-n**2)*dd_pi &
                      -0.5d0*sqrt(3d0)*l*n*(l**2+m**2)*dd_delta !3z2-r2,zx
    E2c_int(9,8,i)  = 0.5d0*sqrt(3d0)*(l**2-m**2)*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                      +sqrt(3d0)*n**2*(m**2-l**2)*dd_pi &
                      +0.25d0*sqrt(3d0)*(1d0+n**2)*(l**2-m**2)*dd_delta !3z2-r2,x2-y2
    E2c_int(9,9,i)  = (n**2-0.5d0*(l**2+m**2))**2*dd_sigma &
                      +3d0*n**2*(l**2+m**2)*dd_pi &
                      +0.75d0*(l**2+m**2)**2*dd_delta !3z2-r2,3z2-r2
    E2c_int(9,10,i) = (n**2-0.5d0*(l**2+m**2))*s2cda_sigma !3z2-r2,s2

! anion s2-orbital
    E2c_int(10,1,i)  = s2a_sc_sigma !s2,s
    E2c_int(10,2,i)  = l*s2apc_sigma !s2,x
    E2c_int(10,3,i)  = m*s2apc_sigma !s2,y
    E2c_int(10,4,i)  = n*s2apc_sigma !s2,z
    E2c_int(10,5,i)  = sqrt(3d0)*l*m*s2adc_sigma !s2,xy
    E2c_int(10,6,i)  = sqrt(3d0)*m*n*s2adc_sigma !s2,yz
    E2c_int(10,7,i)  = sqrt(3d0)*n*l*s2adc_sigma !s2,zx
    E2c_int(10,8,i)  = 0.5d0*sqrt(3d0)*(l**2-m**2)*s2adc_sigma !s2,x2-y2
    E2c_int(10,9,i)  = (n**2-0.5d0*(l**2+m**2))*s2adc_sigma !s2,3z2-r2
    E2c_int(10,10,i) = s2s2_sigma !s2,s2

  end do

end subroutine calc_two_center_integral
!----------------------------------------------------------------------------
subroutine calc_zham_mat
  implicit  none
  integer :: ik,i,ia,ic
  complex(8) :: zphase


  zham_mat = 0d0
  do ik = nk_s, nk_e

    zham_mat( 1, 1,ik)=Ea_s
    zham_mat( 2, 2,ik)=Ea_p
    zham_mat( 3, 3,ik)=Ea_p
    zham_mat( 4, 4,ik)=Ea_p
    zham_mat( 5, 5,ik)=Ea_d
    zham_mat( 6, 6,ik)=Ea_d
    zham_mat( 7, 7,ik)=Ea_d
    zham_mat( 8, 8,ik)=Ea_d
    zham_mat( 9, 9,ik)=Ea_d
    zham_mat(10,10,ik)=Ea_s2

    zham_mat(11,11,ik)=Ec_s
    zham_mat(12,12,ik)=Ec_p
    zham_mat(13,13,ik)=Ec_p
    zham_mat(14,14,ik)=Ec_p
    zham_mat(15,15,ik)=Ec_d
    zham_mat(16,16,ik)=Ec_d
    zham_mat(17,17,ik)=Ec_d
    zham_mat(18,18,ik)=Ec_d
    zham_mat(19,19,ik)=Ec_d
    zham_mat(20,20,ik)=Ec_s2
    do ia = 1, 10
      do ic = 1, 10

        
        do i = 1, num_nearest_neighbor
          zphase=exp(zi*sum(kvec(:,ik)*Rvec_ac(:,i)))
          zham_mat(ia,10+ic,ik)=zham_mat(ia,10+ic,ik) + E2c_int(ia,ic,i)*zphase
        end do

      end do
    end do


    do ia = 1, 10
      do ic = 1, 10
        zham_mat(10+ic,ia,ik)=conjg(zham_mat(ia,10+ic,ik))
      end do
    end do
  end do

  zham_mat(21:40,21:40,:) = zham_mat(1:20,1:20,:)

  zham_mat(2,20+4,:) = zham_mat(2,20+4,:)    + delta_a/3d0  ! <x;up   |SO| z; down>
  zham_mat(20+2,4,:) = zham_mat(20+2,4,:)    - delta_a/3d0  ! <x;down |SO| z; up>
  zham_mat(3,20+4,:) = zham_mat(3,20+4,:) - zi*delta_a/3d0  ! <y;up   |SO| z; down>
  zham_mat(20+3,4,:) = zham_mat(20+3,4,:) - zi*delta_a/3d0  ! <y;down |SO| z; up>
  zham_mat(2,3,:)    = zham_mat(2,3,:)    - zi*delta_a/3d0  ! <x;up   |SO| y; up> !! check
  zham_mat(2+20,3+20,:) = zham_mat(2+20,3+20,:)+zi*delta_a/3d0 ! <x;down   |SO| y; down> !! check

  zham_mat(20+4,2,:) = conjg(zham_mat(2,20+4,:))            ! <z; down|SO| x;up>
  zham_mat(4,20+2,:) = conjg(zham_mat(20+2,4,:))            ! <z; up  |SO| x;down>
  zham_mat(20+4,3,:) = conjg(zham_mat(3,20+4,:))            ! <z; down|SO| y;up>
  zham_mat(4,20+3,:) = conjg(zham_mat(20+3,4,:))            ! <z; up  |SO| y;down>
  zham_mat(3,2,:)    = conjg(zham_mat(2,3,:))               ! <y; up  |SO| x;up> !! check
  zham_mat(3+20,2+20,:) = conjg(zham_mat(2+20,3+20,:))      ! <y; down|SO| x;down> !! check


  zham_mat(12,20+14,:) = zham_mat(12,20+14,:) +delta_c/3d0     ! <x;up   |SO| z; down>
  zham_mat(20+12,14,:) = zham_mat(20+12,14,:) -delta_c/3d0     ! <x;down |SO| z; up>
  zham_mat(13,20+14,:) = zham_mat(13,20+14,:) -zi*delta_c/3d0  ! <y;up   |SO| z; down>
  zham_mat(20+13,14,:) = zham_mat(20+13,14,:) -zi*delta_c/3d0  ! <y;down |SO| z; up>
  zham_mat(12,13,:)    = zham_mat(12,13,:)    - zi*delta_c/3d0 ! <x;up   |SO| y; up> !! check
  zham_mat(12+20,13+20,:) = zham_mat(12+20,13+20,:)+zi*delta_c/3d0 ! <x;down   |SO| y; down> !! check

  zham_mat(20+14,12,:) = conjg(zham_mat(12,20+14,:))           ! <z; down |SO| x;up>
  zham_mat(14,20+12,:) = conjg(zham_mat(20+12,14,:))           ! <z; up |SO| x;down>
  zham_mat(20+14,13,:) = conjg(zham_mat(13,20+14,:))           ! <z; down|SO| y;up>
  zham_mat(14,20+13,:) = conjg(zham_mat(20+13,14,:))           ! <z; up |SO| y;down>
  zham_mat(13,12,:)    = conjg(zham_mat(12,13,:))               ! <y; up  |SO| x;up> !! check
  zham_mat(13+20,12+20,:) = conjg(zham_mat(12+20,13+20,:))      ! <y; down|SO| x;down> !! check

end subroutine calc_zham_mat
!----------------------------------------------------------------------------
subroutine calc_zJ_mat
  implicit none
  integer :: ik,i,ia,ic
  complex(8) :: zphase

  zJ_mat = 0d0

  do ik = nk_s, nk_e

    do ia = 1, 10
      do ic = 1, 10
        do i = 1, num_nearest_neighbor
          zphase=exp(zi*sum(kvec(:,ik)*Rvec_ac(:,i)))
          zJ_mat(ia,10+ic,1,ik)=zJ_mat(ia,10+ic,1,ik) &
              + E2c_int(ia,ic,i)*zphase*zi*Rvec_ac(1,i)

          zJ_mat(ia,10+ic,2,ik)=zJ_mat(ia,10+ic,2,ik) &
              + E2c_int(ia,ic,i)*zphase*zi*Rvec_ac(2,i)

          zJ_mat(ia,10+ic,3,ik)=zJ_mat(ia,10+ic,3,ik) &
              + E2c_int(ia,ic,i)*zphase*zi*Rvec_ac(3,i)

        end do
      end do
    end do

    do ia = 1, 10
      do ic = 1, 10
        zJ_mat(10+ic,ia,:,ik)=conjg(zJ_mat(ia,10+ic,:,ik))
      end do
    end do

  end do

  zJ_mat(21:40,21:40,:,:) = zJ_mat(1:20,1:20,:,:)

end subroutine calc_zJ_mat
!----------------------------------------------------------------------------
subroutine calc_current(Act_t, jt_t)
  implicit none
  real(8),intent(in) :: Act_t(3)
  real(8),intent(out) :: jt_t(3)
  real(8),allocatable :: Amat_tmp(:,:,:)
  integer :: ik, ib

  jt_t = 0d0
  

  do ik = nk_s, nk_e
    kvec(:,ik) = kvec0(:,ik) + alpha_fsa(ik)**(1d0/npower_focal_spot_average)*Act_t(:)
  end do

  call calc_two_center_integral
  call calc_zJ_mat


  allocate(Amat_tmp(nband,nband,3))
  do ik = nk_s, nk_e

    Amat_tmp(:,:,1) = matmul(zJ_mat(:,:,1,ik),zrho_dm(:,:,ik))
    Amat_tmp(:,:,2) = matmul(zJ_mat(:,:,2,ik),zrho_dm(:,:,ik))
    Amat_tmp(:,:,3) = matmul(zJ_mat(:,:,3,ik),zrho_dm(:,:,ik))

    do ib = 1, nband
      jt_t(1) = jt_t(1) + Amat_tmp(ib,ib,1)/(npower_focal_spot_average*alpha_fsa(ik))
      jt_t(2) = jt_t(2) + Amat_tmp(ib,ib,2)/(npower_focal_spot_average*alpha_fsa(ik))
      jt_t(3) = jt_t(3) + Amat_tmp(ib,ib,3)/(npower_focal_spot_average*alpha_fsa(ik))
    end do

  end do

  call comm_allreduce(jt_t)
  jt_t = jt_t/(nkpoint*volume)


end subroutine calc_current
!----------------------------------------------------------------------------
subroutine calc_bandstructure_zincblende
  implicit none
  integer :: ik, nkpoint_8
  real(8):: kvec_L(3), kvec_Gamma(3), kvec_X(3)
  real(8):: kvec_W(3), kvec_K(3), kvec_U(3)
  real(8),allocatable :: kk(:)
!LAPACK
  integer :: ndim
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

  ndim = 40
  lwork = 4*ndim**2+4*ndim+256
  allocate(work_lp(lwork),rwork(3*ndim-2),w(ndim))


!LAPACK
  allocate(kk(nkpoint)); kk = 0d0

  if(.not.if_root_global)return

  deallocate(zham_mat)
  allocate(zham_mat(1:ndim,1:ndim,1:nkpoint))

!  kvec_L=2d0*pi/lattice_const*0.5d0
!  kvec_Gamma = 0d0
!  kvec_X=0d0; kvec_X(1)=2d0*pi/lattice_const
!  kvec_W=0d0; kvec_W(1)=2d0*pi/lattice_const; kvec_W(2)=2d0*pi/lattice_const*0.5d0
!  kvec_K=2d0*pi/lattice_const*0.75d0; kvec_K(3)=0d0
!  kvec_U(1)=1d0;kvec_U(2)=0.25d0;kvec_U(3)=0.25d0; kvec_U=kvec_U*2d0*pi/lattice_const

  kvec_Gamma = 0d0
  kvec_L = 0.50d0*reciprocal_lattice_vec(:,1) &
          +0.50d0*reciprocal_lattice_vec(:,2) &
          +0.50d0*reciprocal_lattice_vec(:,3)

  kvec_X = 0.00d0*reciprocal_lattice_vec(:,1) &
          +0.50d0*reciprocal_lattice_vec(:,2) &
          +0.50d0*reciprocal_lattice_vec(:,3)

  kvec_W = 0.25d0*reciprocal_lattice_vec(:,1) &
          +(3d0/4d0)*reciprocal_lattice_vec(:,2) &
          +0.50d0*reciprocal_lattice_vec(:,3)

  kvec_K = (3d0/8d0)*reciprocal_lattice_vec(:,1) &
          +(3d0/4d0)*reciprocal_lattice_vec(:,2) &
          +(3d0/8d0)*reciprocal_lattice_vec(:,3)

  kvec_U = (1d0/4d0)*reciprocal_lattice_vec(:,1) &
          +(5d0/8d0)*reciprocal_lattice_vec(:,2) &
          +(5d0/8d0)*reciprocal_lattice_vec(:,3)


  kvec0 = 0d0
  
  nkpoint_8 = nkpoint/9
  if(mod(nkpoint,9)/=0)then
    write(*,*)"mod(nkpoint,8)/=0"
    stop
  end if

!L to Gamma
  do ik = 1,nkpoint_8
    kvec(:,ik+0*nkpoint_8) = kvec_L(:) +(ik-1)*(kvec_Gamma(:)-kvec_L(:))/nkpoint_8
    if(ik /= 1)kk(ik+0*nkpoint_8)=kk(ik+0*nkpoint_8-1) &
        +sqrt(sum((kvec_Gamma(:)-kvec_L(:))**2))/nkpoint_8
  end do

!Gamma to X
  do ik = 1,nkpoint_8
    kvec(:,ik+1*nkpoint_8) = kvec_Gamma(:) +(ik-1)*(kvec_X(:)-kvec_Gamma(:))/nkpoint_8
    if(ik==1)then
      kk(ik+1*nkpoint_8)=kk(ik+1*nkpoint_8-1)+sqrt(sum((kvec_Gamma(:)-kvec_L(:))**2))/nkpoint_8
    else
      kk(ik+1*nkpoint_8)=kk(ik+1*nkpoint_8-1)+sqrt(sum((kvec_X(:)-kvec_Gamma(:))**2))/nkpoint_8
    end if
  end do

!X to W
  do ik = 1,nkpoint_8
    kvec(:,ik+2*nkpoint_8) = kvec_X(:) +(ik-1)*(kvec_W(:)-kvec_X(:))/nkpoint_8
    if(ik==1)then
      kk(ik+2*nkpoint_8)=kk(ik+2*nkpoint_8-1)+sqrt(sum((kvec_X(:)-kvec_Gamma(:))**2))/nkpoint_8
    else
      kk(ik+2*nkpoint_8)=kk(ik+2*nkpoint_8-1)+sqrt(sum((kvec_W(:)-kvec_X(:))**2))/nkpoint_8
    end if
  end do

!W to K
  do ik = 1,nkpoint_8
    kvec(:,ik+3*nkpoint_8) = kvec_W(:) +(ik-1)*(kvec_K(:)-kvec_W(:))/nkpoint_8
    if(ik==1)then
      kk(ik+3*nkpoint_8)=kk(ik+3*nkpoint_8-1)+sqrt(sum((kvec_W(:)-kvec_X(:))**2))/nkpoint_8
    else
      kk(ik+3*nkpoint_8)=kk(ik+3*nkpoint_8-1)+sqrt(sum((kvec_K(:)-kvec_W(:))**2))/nkpoint_8
    end if
  end do

!K to L
  do ik = 1,nkpoint_8
    kvec(:,ik+4*nkpoint_8) = kvec_K(:) +(ik-1)*(kvec_L(:)-kvec_K(:))/nkpoint_8
    if(ik==1)then
      kk(ik+4*nkpoint_8)=kk(ik+4*nkpoint_8-1)+sqrt(sum((kvec_K(:)-kvec_W(:))**2))/nkpoint_8
    else
      kk(ik+4*nkpoint_8)=kk(ik+4*nkpoint_8-1)+sqrt(sum((kvec_L(:)-kvec_K(:))**2))/nkpoint_8
    end if
  end do

!L to W
  do ik = 1,nkpoint_8
    kvec(:,ik+5*nkpoint_8) = kvec_L(:) +(ik-1)*(kvec_W(:)-kvec_L(:))/nkpoint_8
    if(ik==1)then
      kk(ik+5*nkpoint_8)=kk(ik+5*nkpoint_8-1)+sqrt(sum((kvec_L(:)-kvec_K(:))**2))/nkpoint_8
    else
      kk(ik+5*nkpoint_8)=kk(ik+5*nkpoint_8-1)+sqrt(sum((kvec_W(:)-kvec_L(:))**2))/nkpoint_8
    end if
  end do

!W to X
  do ik = 1,nkpoint_8
    kvec(:,ik+6*nkpoint_8) = kvec_W(:) +(ik-1)*(kvec_X(:)-kvec_W(:))/nkpoint_8
    if(ik==1)then
      kk(ik+6*nkpoint_8)=kk(ik+6*nkpoint_8-1)+sqrt(sum((kvec_W(:)-kvec_L(:))**2))/nkpoint_8
    else
      kk(ik+6*nkpoint_8)=kk(ik+6*nkpoint_8-1)+sqrt(sum((kvec_X(:)-kvec_W(:))**2))/nkpoint_8
    end if
  end do

!X to U
  do ik = 1,nkpoint_8
    kvec(:,ik+7*nkpoint_8) = kvec_X(:) +(ik-1)*(kvec_U(:)-kvec_X(:))/nkpoint_8
    if(ik==1)then
      kk(ik+7*nkpoint_8)=kk(ik+7*nkpoint_8-1)+sqrt(sum((kvec_X(:)-kvec_W(:))**2))/nkpoint_8
    else
      kk(ik+7*nkpoint_8)=kk(ik+7*nkpoint_8-1)+sqrt(sum((kvec_U(:)-kvec_X(:))**2))/nkpoint_8
    end if
  end do


!K to Gamma
  do ik = 1,nkpoint_8
    kvec(:,ik+8*nkpoint_8) = kvec_K(:) +(ik-1)*(kvec_Gamma(:)-kvec_K(:))/(nkpoint_8-1)
    if(ik==1)then
      kk(ik+8*nkpoint_8)=kk(ik+8*nkpoint_8-1)+sqrt(sum((kvec_U(:)-kvec_X(:))**2))/nkpoint_8
    else
      kk(ik+8*nkpoint_8)=kk(ik+8*nkpoint_8-1)+sqrt(sum((kvec_Gamma(:)-kvec_K(:))**2))/(nkpoint_8-1)
    end if
  end do


!  kvec(:,1)=kvec_K
!  kvec(:,2)=kvec_U

  call calc_two_center_integral
  call calc_zham_mat


  open(20,file="bandstructure.out")
  do ik = 1, nkpoint

    call zheev('V', 'U', ndim, zham_mat(1:ndim,1:ndim,ik), ndim, w, work_lp, lwork, rwork, info)
    write(20,"(999e26.16e3)")kk(ik)/kk(nkpoint),kvec(:,ik),w

  end do
  close(20)

  deallocate(zham_mat)
  allocate(zham_mat(1:ndim,1:ndim,nk_s:nk_e))

end subroutine calc_bandstructure_zincblende
!----------------------------------------------------------------------------
subroutine set_equilibrium_density_matrix
  implicit none
  integer :: ik
  complex(8),allocatable :: zAmat_tmp(:,:)
!LAPACK
  integer :: ndim
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

  ndim = 40
  lwork = 4*ndim**2+4*ndim+256
  allocate(work_lp(lwork),rwork(3*ndim-2),w(ndim))


!LAPACK


  allocate(zAmat_tmp(nband,nband))

  if(.not. allocated(zrho_dm))then
    allocate(zrho_dm(nband,nband,nk_s:nk_e))
  end if

  kvec = kvec0
  call calc_two_center_integral
  call calc_zham_mat

  zrho_dm = 0d0
  do ik = nk_s, nk_e
    call zheev('V', 'U', ndim, zham_mat(1:ndim,1:ndim,ik), ndim, w, work_lp, lwork, rwork, info)

    zrho_dm(1,1,ik)=1d0; zrho_dm(2,2,ik)=1d0; zrho_dm(3,3,ik)=1d0; zrho_dm(4,4,ik)=1d0
    zrho_dm(5,5,ik)=1d0; zrho_dm(6,6,ik)=1d0; zrho_dm(7,7,ik)=1d0; zrho_dm(8,8,ik)=1d0

    zAmat_tmp = matmul(matmul(zham_mat(:,:,ik),zrho_dm(:,:,ik)), &
                       conjg(transpose(zham_mat(:,:,ik))))
    zrho_dm(:,:,ik) = zAmat_tmp

  end do

end subroutine set_equilibrium_density_matrix
!----------------------------------------------------------------------------
subroutine dt_evolve_elec_system(Act_in,dt_in)
  implicit none
  real(8),intent(in) :: Act_in(3), dt_in
  integer :: ik, ib, ib1,ib2
  complex(8),allocatable :: zAmat_tmp(:,:),zBmat_tmp(:,:)
  real(8),allocatable :: ref_pop(:)
!LAPACK
  integer :: ndim
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

  ndim = 40
  lwork = 4*ndim**2+4*ndim+256
  allocate(work_lp(lwork),rwork(3*ndim-2),w(ndim))


!LAPACK


  do ik = nk_s, nk_e
    kvec(:,ik) = kvec0(:,ik) + alpha_fsa(ik)**(1d0/npower_focal_spot_average)*Act_in(:)
  end do
  allocate(zAmat_tmp(nband,nband),zBmat_tmp(nband,nband))
  allocate(ref_pop(nband))
  ref_pop(1:8) = 1d0; ref_pop(9:nband) = 0d0

  call calc_two_center_integral
  call calc_zham_mat

  do ik = nk_s, nk_e
    call zheev('V', 'U', ndim, zham_mat(1:ndim,1:ndim,ik), ndim, w, work_lp, lwork, rwork, info)
    zAmat_tmp = matmul( &
        matmul(conjg(transpose(zham_mat(:,:,ik))),zrho_dm(:,:,ik)) &
        ,zham_mat(:,:,ik))

    do ib = 1, nband
      zAmat_tmp(ib,ib)=zAmat_tmp(ib,ib)-ref_pop(ib)
    end do
    
    do ib1 = 1, nband
      zBmat_tmp(ib1,ib1) = exp(-dt_in/T1_relax)
      do ib2 = ib1+1,nband
        zBmat_tmp(ib1,ib2) = exp(-zi*(w(ib1)-w(ib2))*dt_in -dt_in/T2_relax)
        zBmat_tmp(ib2,ib1) = exp(-zi*(w(ib2)-w(ib1))*dt_in -dt_in/T2_relax)
      end do
    end do

    zAmat_tmp = zAmat_tmp*zBmat_tmp 
    do ib = 1, nband
      zAmat_tmp(ib,ib)=zAmat_tmp(ib,ib)+ref_pop(ib)
    end do

    zrho_dm(:,:,ik) = matmul( &
        matmul(zham_mat(:,:,ik),zAmat_tmp) &
        ,conjg(transpose(zham_mat(:,:,ik))))


  end do

end subroutine dt_evolve_elec_system
!----------------------------------------------------------------------------
subroutine dt_evolve_elec_system_mod(Act_1_in, Act_2_in, dt_in)
  implicit none
  real(8),intent(in) :: Act_1_in(3), Act_2_in(3), dt_in
  integer :: ik, ib, ib1,ib2
  complex(8),allocatable :: zAmat_tmp(:,:),zBmat_tmp(:,:)
  real(8),allocatable :: ref_pop(:)
  real(8),allocatable :: eps_1(:),eps_2(:)
  complex(8),allocatable :: zham_mat_1(:,:,:),zham_mat_2(:,:,:)
  complex(8),allocatable :: zUm(:,:)
!LAPACK
  integer :: ndim
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

  ndim = 40
  lwork = 4*ndim**2+4*ndim+256
  allocate(work_lp(lwork),rwork(3*ndim-2),w(ndim))
 

!LAPACK
  allocate(eps_1(ndim),eps_2(ndim))
  allocate(zham_mat_1(nband,nband,nk_s:nk_e),zham_mat_2(nband,nband,nk_s:nk_e))
  allocate(zAmat_tmp(nband,nband),zBmat_tmp(nband,nband))
  allocate(ref_pop(nband))
  ref_pop(1:8) = 1d0; ref_pop(9:nband) = 0d0
  allocate(zUm(ndim,ndim))


  do ik = nk_s, nk_e
    kvec(:,ik) = kvec0(:,ik) + alpha_fsa(ik)**(1d0/npower_focal_spot_average)*Act_1_in(:)
  end do
  call calc_two_center_integral
  call calc_zham_mat
  zham_mat_1 = zham_mat

  do ik = nk_s, nk_e
    kvec(:,ik) = kvec0(:,ik) + alpha_fsa(ik)**(1d0/npower_focal_spot_average)*Act_2_in(:)
  end do
  call calc_two_center_integral
  call calc_zham_mat
  zham_mat_2 = zham_mat


  do ik = nk_s, nk_e
    call zheev('V', 'U', ndim, zham_mat_1(1:ndim,1:ndim,ik), ndim, eps_1, work_lp, lwork, rwork, info)
    call zheev('V', 'U', ndim, zham_mat_2(1:ndim,1:ndim,ik), ndim, eps_2, work_lp, lwork, rwork, info)

    do ib1 = 1,ndim
      do ib2 = 1,ndim
        zUm(ib2,ib1)=sum(conjg(zham_mat_2(:,ib2,ik))*zham_mat_1(:,ib1,ik)) &
            *exp(-0.5d0*zi*dt_in*(eps_2(ib2)+eps_1(ib1)))
      end do
    end do
    
! blocking to freeze the bands
    zUm = zUm*blocking_matrix_band_frozen

! convert to the H1 basis expression
    zAmat_tmp = matmul( &
        matmul(conjg(transpose(zham_mat_1(:,:,ik))),zrho_dm(:,:,ik)) &
        ,zham_mat_1(:,:,ik))

! apply the first relaxation at t
    do ib = 1, nband
      zAmat_tmp(ib,ib)=zAmat_tmp(ib,ib)-ref_pop(ib)
    end do
    
    do ib1 = 1, nband
      zBmat_tmp(ib1,ib1) = exp(-0.5d0*dt_in/T1_relax)
      do ib2 = ib1+1,nband
        zBmat_tmp(ib1,ib2) = exp(-0.5d0*dt_in/T2_relax)
        zBmat_tmp(ib2,ib1) = exp(-0.5d0*dt_in/T2_relax)
      end do
    end do

    zBmat_tmp = zAmat_tmp*zBmat_tmp 
    do ib = 1, nband
      zBmat_tmp(ib,ib)=zBmat_tmp(ib,ib)+ref_pop(ib)
    end do


! apply the time propagation from t to t+dt
    zAmat_tmp = matmul(zUm(:,:),matmul(zBmat_tmp, conjg(transpose(zUm(:,:)))))

! apply the second relaxation at t+dt
    do ib = 1, nband
      zAmat_tmp(ib,ib)=zAmat_tmp(ib,ib)-ref_pop(ib)
    end do
    
    do ib1 = 1, nband
      zBmat_tmp(ib1,ib1) = exp(-0.5d0*dt_in/T1_relax)
      do ib2 = ib1+1,nband
        zBmat_tmp(ib1,ib2) = exp(-0.5d0*dt_in/T2_relax)
        zBmat_tmp(ib2,ib1) = exp(-0.5d0*dt_in/T2_relax)
      end do
    end do

    zBmat_tmp = zAmat_tmp*zBmat_tmp 
    do ib = 1, nband
      zBmat_tmp(ib,ib)=zBmat_tmp(ib,ib)+ref_pop(ib)
    end do

! conver to the original basis expression
    zrho_dm(:,:,ik) = matmul( &
        matmul(zham_mat_2(:,:,ik),zBmat_tmp) &
        ,conjg(transpose(zham_mat_2(:,:,ik))))


  end do

end subroutine dt_evolve_elec_system_mod
!----------------------------------------------------------------------------
subroutine calc_num_electron(num_elec)
  implicit none
  real(8),intent(out) :: num_elec
  integer :: ik, ib


  num_elec = 0d0
  do ik = nk_s, nk_e
    do ib = 1, nband
      num_elec = num_elec + zrho_dm(ib,ib,ik)
    end do
  end do

  call comm_allreduce(num_elec)
  num_elec = num_elec/(nkpoint)
  
end subroutine calc_num_electron
!----------------------------------------------------------------------------
subroutine prepare_sampling_for_focal_spot_average
  implicit none
  integer :: isample
  integer :: ik
  real(8) :: x1,x2,x3,alpha

  ik = 0
  do isample = 1, num_sample_focal_spot_average

    call vdCorput_sequence(isample,2,x1); x1 = x1-0.5d0
    call vdCorput_sequence(isample,3,x2); x2 = x2-0.5d0
    call vdCorput_sequence(isample,5,x3); x3 = x3-0.5d0
    call vdCorput_sequence(isample,7,alpha)

! x1,x2,x3,+
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =( x1)*reciprocal_lattice_vec(:,1) &
                  +( x2)*reciprocal_lattice_vec(:,2) &
                  +( x3)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if

! x1,x2,x3,-
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =(-x1)*reciprocal_lattice_vec(:,1) &
                  +(-x2)*reciprocal_lattice_vec(:,2) &
                  +(-x3)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if

! x1,x3,x2,+
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =( x1)*reciprocal_lattice_vec(:,1) &
                  +( x3)*reciprocal_lattice_vec(:,2) &
                  +( x2)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if

! x1,x3,x2,-
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =(-x1)*reciprocal_lattice_vec(:,1) &
                  +(-x3)*reciprocal_lattice_vec(:,2) &
                  +(-x2)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if


! x2,x1,x3,+
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =( x2)*reciprocal_lattice_vec(:,1) &
                  +( x1)*reciprocal_lattice_vec(:,2) &
                  +( x3)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if

! x2,x1,x3,-
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =(-x2)*reciprocal_lattice_vec(:,1) &
                  +(-x1)*reciprocal_lattice_vec(:,2) &
                  +(-x3)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if

! x2,x3,x1,+
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =( x2)*reciprocal_lattice_vec(:,1) &
                  +( x3)*reciprocal_lattice_vec(:,2) &
                  +( x1)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if

! x2,x3,x1,-
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =(-x2)*reciprocal_lattice_vec(:,1) &
                  +(-x3)*reciprocal_lattice_vec(:,2) &
                  +(-x1)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if

! x3,x1,x2,+
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =( x3)*reciprocal_lattice_vec(:,1) &
                  +( x1)*reciprocal_lattice_vec(:,2) &
                  +( x2)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if

! x3,x1,x2,-
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =(-x3)*reciprocal_lattice_vec(:,1) &
                  +(-x1)*reciprocal_lattice_vec(:,2) &
                  +(-x2)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if

! x3,x2,x1,+
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =( x3)*reciprocal_lattice_vec(:,1) &
                  +( x2)*reciprocal_lattice_vec(:,2) &
                  +( x1)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if

! x3,x2,x1,-
    ik = ik + 1
    if(ik >= nk_s .and. ik <= nk_e)then
      kvec0(:,ik) =(-x3)*reciprocal_lattice_vec(:,1) &
                  +(-x2)*reciprocal_lattice_vec(:,2) &
                  +(-x1)*reciprocal_lattice_vec(:,3) 

      alpha_fsa(ik) = alpha
    end if


  end do

end subroutine prepare_sampling_for_focal_spot_average
!----------------------------------------------------------------------------
subroutine set_blocking_matrix_band_frozen
  implicit none
  integer :: ib1, ib2
  integer,allocatable :: nflag_include(:)

  allocate(nflag_include(nband))

  nflag_include( 1) = 1; nflag_include( 2) = 1
  nflag_include( 3) = 1; nflag_include( 4) = 1
  nflag_include( 5) = 1; nflag_include( 6) = 1
  nflag_include( 7) = 1; nflag_include( 8) = 1  ! valence top

  nflag_include( 9) = 1; nflag_include(10) = 1  ! conduction bottom
  nflag_include(11) = 1; nflag_include(12) = 1
  nflag_include(13) = 1; nflag_include(14) = 1
  nflag_include(15) = 1; nflag_include(16) = 1
  nflag_include(17) = 1; nflag_include(18) = 1
  nflag_include(19) = 1; nflag_include(20) = 1
  nflag_include(21) = 1; nflag_include(22) = 1
  nflag_include(23) = 1; nflag_include(24) = 1
  nflag_include(25) = 1; nflag_include(26) = 1
  nflag_include(27) = 1; nflag_include(28) = 1
  nflag_include(29) = 1; nflag_include(30) = 1
  nflag_include(31) = 1; nflag_include(32) = 1
  nflag_include(33) = 1; nflag_include(34) = 1
  nflag_include(35) = 1; nflag_include(36) = 1
  nflag_include(37) = 1; nflag_include(38) = 1
  nflag_include(39) = 1; nflag_include(40) = 1


  blocking_matrix_band_frozen = 1d0
! include only the bottom of conduction and top of valence bands
  do ib1 = 1, nband
    do ib2 = 1, nband
      if(ib1 /= ib2)then
        blocking_matrix_band_frozen(ib1,ib2)=nflag_include(ib1)*nflag_include(ib2)
      end if
    end do
  end do

end subroutine set_blocking_matrix_band_frozen
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

end module electronic_system
