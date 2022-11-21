module global_variables
  use parallel
  use communication
  use math
  use constants
  use inputoutput

! material class
  integer,parameter :: n_diamond = 0, n_zincblende = 1
  integer,parameter :: material_class = n_zincblende

! tight-binding
  real(8) :: lattice_vec(3,3), lattice_const, volume
  real(8) :: reciprocal_lattice_vec(3,3)
  integer,parameter :: nspin = 2
  integer :: nband, nkpoint
  integer :: nkx, nky, nkz
  real(8),allocatable :: kvec0(:,:),kvec(:,:)
  complex(8),allocatable :: zpsi(:,:,:)
  complex(8),allocatable :: zHam_mat(:,:,:)

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



end module global_variables
