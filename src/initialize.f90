subroutine initialize
  use global_variables
  implicit none


  allocate(zpsi(nband,nspin,nkpoint), zHam_mat(nband,nband,nkpoint))

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
!! spin down components (31-40)

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


  num_nearest_neighbor = 4
  allocate(Rvec_ac(3,num_nearest_neighbor))
  allocate(E2c_int(10,10,num_nearest_neighbor))

  Rvec_ac(:,1) = lattice_const*0.25d0
  Rvec_ac(:,2) = lattice_const*0.25d0 - lattice_vec(:,1)
  Rvec_ac(:,3) = lattice_const*0.25d0 - lattice_vec(:,2)
  Rvec_ac(:,4) = lattice_const*0.25d0 - lattice_vec(:,3)
! check distance
  write(*,*)"dist 1=",sqrt(sum(Rvec_ac(:,1)**2))/lattice_const
  write(*,*)"dist 2=",sqrt(sum(Rvec_ac(:,2)**2))/lattice_const
  write(*,*)"dist 3=",sqrt(sum(Rvec_ac(:,3)**2))/lattice_const
  write(*,*)"dist 4=",sqrt(sum(Rvec_ac(:,4)**2))/lattice_const


  include "include_tb_parameters/set_GaAs_Jancu1998.f90"

end subroutine initialize
