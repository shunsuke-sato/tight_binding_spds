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
    E2c_int(2,1,i) = (-1d0)*l*scpa_sigma !xs

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
    E2c_int(3,4,i) = m*n*pp_sigma-l*m*pp_pi !y,z
    E2c_int(3,5,i) = sqrt(3d0)*m**2*l*padc_sigma+l*(1d0-2d0*m**2)*padc_pi !y,xy
    E2c_int(3,6,i) = sqrt(3d0)*m**2*n*padc_sigma+n*(1d0-2d0*m**2)*padc_pi !y,yz
    E2c_int(3,7,i) = sqrt(3d0)*m*n*l*padc_sigma-2d0*m*n*l*padc_pi !y,zx
    E2c_int(3,8,i) = 0.5d0*sqrt(3d0)*m*(l**2-m**2)*padc_sigma &
                     -m*(1d0+l**2-m**2)*padc_pi !y,x2-y2
    E2c_int(3,9,i) = m*(n**2-0.5d0*(l**2+m**2))*padc_sigma &
                     -sqrt(3d0)*m*n**2*padc*pi !y,3r2-z2
    E2c_int(3,10,i) = (-1d0)*m*s2cpa_sigma !y,s2


  end do

end subroutine calc_two_center_integral
