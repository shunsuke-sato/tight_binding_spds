subroutine initialize
  use global_variables
  implicit none
  real(8) :: vec_t(3)

  allocate(kvec0(3,nkpoint),kvec(3,nkpoint))
  kvec0 = 0d0; kvec=0d0
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
  volume = sum(lattice_vec(:,1)*vec_t(:))
  write(*,*)"volume=",volume

  reciprocal_lattice_vec(:,1)=2d0*pi/volume*cross_product(lattice_vec(:,2),lattice_vec(:,3))
  reciprocal_lattice_vec(:,2)=2d0*pi/volume*cross_product(lattice_vec(:,3),lattice_vec(:,1))
  reciprocal_lattice_vec(:,3)=2d0*pi/volume*cross_product(lattice_vec(:,1),lattice_vec(:,2))

  write(*,*)"a1*b123/2pi",sum(lattice_vec(:,1)*reciprocal_lattice_vec(:,1))/(2d0*pi) &
                         ,sum(lattice_vec(:,1)*reciprocal_lattice_vec(:,2))/(2d0*pi) &
                         ,sum(lattice_vec(:,1)*reciprocal_lattice_vec(:,3))/(2d0*pi)

  write(*,*)"a2*b123/2pi",sum(lattice_vec(:,2)*reciprocal_lattice_vec(:,1))/(2d0*pi) &
                         ,sum(lattice_vec(:,2)*reciprocal_lattice_vec(:,2))/(2d0*pi) &
                         ,sum(lattice_vec(:,2)*reciprocal_lattice_vec(:,3))/(2d0*pi)

  write(*,*)"a3*b123/2pi",sum(lattice_vec(:,3)*reciprocal_lattice_vec(:,1))/(2d0*pi) &
                         ,sum(lattice_vec(:,3)*reciprocal_lattice_vec(:,2))/(2d0*pi) &
                         ,sum(lattice_vec(:,3)*reciprocal_lattice_vec(:,3))/(2d0*pi)

  num_nearest_neighbor = 4
  allocate(Rvec_ac(3,num_nearest_neighbor))
  allocate(E2c_int(10,10,num_nearest_neighbor))

!  Rvec_ac(:,1) = lattice_const*0.25d0
  Rvec_ac(:,1) = lattice_const*0.25d0
  Rvec_ac(:,2) = Rvec_ac(:,1) - lattice_vec(:,1)
  Rvec_ac(:,3) = Rvec_ac(:,1) - lattice_vec(:,2)
  Rvec_ac(:,4) = Rvec_ac(:,1) - lattice_vec(:,3)
! check distance
  write(*,*)"dist 1=",sqrt(sum(Rvec_ac(:,1)**2))/lattice_const
  write(*,*)"dist 2=",sqrt(sum(Rvec_ac(:,2)**2))/lattice_const
  write(*,*)"dist 3=",sqrt(sum(Rvec_ac(:,3)**2))/lattice_const
  write(*,*)"dist 4=",sqrt(sum(Rvec_ac(:,4)**2))/lattice_const


  write(*,*)"vec1=",Rvec_ac(:,1)
  write(*,*)"vec1=",Rvec_ac(:,2)
  write(*,*)"vec1=",Rvec_ac(:,3)
  write(*,*)"vec1=",Rvec_ac(:,4)

!  include "include_tb_parameters/set_GaAs_Jancu1998.f90"
  include "include_tb_parameters/set_GaAs_Tan2013.f90"

contains
  function cross_product(vec1, vec2)
    real(8) :: cross_product(3)
    real(8) :: vec1(3), vec2(3)

    cross_product(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
    cross_product(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
    cross_product(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)

  end function cross_product

end subroutine initialize
