subroutine calc_two_center_integral
  use global_variables
  implicit none
  real(8) :: l,m,n
  real(8) :: li,mi,ni
  integer :: i

  E2c_int = 0d0

  do i = 1, num_nearest_neighbor
    l = Rvec_ac(1,i)/sqrt(sum(Rvec_ac(:,i)**2))
    m = Rvec_ac(2,i)/sqrt(sum(Rvec_ac(:,i)**2))
    n = Rvec_ac(3,i)/sqrt(sum(Rvec_ac(:,i)**2))
    li = -l
    mi = -m
    ni = -n

    E2c_int(1,1,i) = ss_sigma !ss
    E2c_int(1,2,i) = l*sapc_sigma !sx
    E2c_int(1,3,i) = m*sapc_sigma !sy
    E2c_int(1,4,i) = n*sapc_sigma !sz

    E2c_int(1,5,i) = sqrt(3d0)*l*m*sadc_sigma !s,xy
    E2c_int(1,6,i) = sqrt(3d0)*m*n*sadc_sigma !s,yx
    E2c_int(1,7,i) = sqrt(3d0)*n*l*sadc_sigma !s,zx
    E2c_int(1,8,i) = 0.5d0*sqrt(3d0)*(l**2-m**2)*sadc_sigma !s,x2-y2
    E2c_int(1,9,i) = (n**2-0.5d0*(l**2+m**2))*sadc_sigma !s,3r2-z2

    E2c_int(1,10,i) = sa_s2c_sigma !ss2


  end do

end subroutine calc_two_center_integral
