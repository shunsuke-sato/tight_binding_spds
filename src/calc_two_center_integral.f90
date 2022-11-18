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
                     -sqrt(3d0)*m*n**2*padc*pi !y,3z2-r2
    E2c_int(3,10,i) = (-1d0)*m*s2cpa_sigma !y,s2

! anion pz-orbital
    E2c_int(4,1,i)  = (-1d0)*n*scpa_sigma !z,s
    E2c_int(4,2,i)  = n*l*pp_sigma-n*l*pp_pi !z,x
    E2c_int(4,3,i)  = m*n*pp_sigma-m*n*pp_pi !z,y
    E2c_int(4,4,i)  = n**2*pp_sigma+(1d0-n**2)*p_pi !z,z
    E2c_int(4,5,i)  = sqrt(3d0)*l*m*n*padc_sigma-2d0*l*m*n*padc_pi !z,xy
    E2c_int(4,6,i)  = sqrt(3d0)*n**2*m*padc_sigma+m*(1d0-2d0*n**2)*padc_pi !z,yz
    E2c_int(4,7,i)  = sqrt(3d0)*n**2*l*padc_sigma+l*(1d0-2d0*n**2)*padc_pi !z,zx
    E2c_int(4,8,i)  = 0.5d0*sqrt(3d0)*n*(l**2-m**2)*padc_sigma &
                      -n**(l**2-m**2)*padc_pi !z,x2-y2
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
                       0.5d0*l*m*(l**2-m**2)*dd_delta !xy,x2-y2
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
                      +m*(1d0-1d0*m**2)*pcda_pi) !yz,z
    E2c_int(6,5,i)  = 3d0*l*m**2*n*dd_sigma+l*n*(1d0-4d0*m**2)*dd_pi &
                      l*n*(m**2-1d0)*dd_delta !yz,xy
    E2c_int(6,6,i)  = 3d0*m**2*n**2*dd_sigma+(m**2+n**2-4d0*m**2*n**2)*dd_pi &
                      +(l**2+m**2*n**2)*dd_delta !yz,yz
    E2c_int(6,7,i)  = 3d0*m*n**2*l*dd_sigma+m*l*(1d0-4d0*n**2)*dd_pi &
                      +m*l*(n**2-1d0)*dd_delta !yz,zx
    E2c_int(6,8,i)  = 1.5d0*m*n*(l**2-m**2)*dd_sigma &
                      -m*n*(1d0+2d0*(l**2-m**2))*dd_pi &
                      +m*n*(1d0+0.5d0*(l**2-m**2))*dd_delta !yz,x2-y2
    E2c_int(6,9,i)  = sqrt(3d0)*m*n*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                      +sqrt(3d0)*m*n*(l**2+m**2-n**2)*dd_pi &
                      -0.5d0*sqrt(3d0)*m*n*(l**2+m**2) !yz,3z2-r2
    E2c_int(6,10,i) = sqrt(3d0)*m*n*s2cda_sigma !yz,s2

! anion dzx-orbital
    E2c_int(7,1,i)  = sqrt(3d0)*n*l*scda_sigma  !zx,s
    E2c_int(7,2,i)  = (-1d0)*(sqrt(3d0)*l**2*n*pcda_sigma &
                      +n*(1d0-2d0*l**2)*pcda_pi) !zx,x
    E2c_int(7,3,i)  = (-1d0)*(sqrt(3d0)*l*m*n*pcda_sigma &
                      -2d0*l*m*n*pcda_pi) !zx,y
    E2c_int(7,4,i)  = (-1d0)*(n**2*l*pcda_sigma &
                      +x*(1d0-2d0*n**2)*pcda_pi) !zx,z
    E2c_int(7,5,i)  = 3d0*l**2*m*n*dd_sigma+m*n*(1d0-4d0*l**2)*dd_pi &
                      +m*n*(l**2-1d0)*dd_delta !zx,xy
    E2c_int(7,6,i)  = 3d0*n**2*l*m*dd_sigma+l*m*(1d0-4d0*n**2)*dd_pi &
                      +l*m*(n**2-1d0)}dd_delta !zx,yz
    E2c_int(7,7,i)  = 3d0*n**2*l**2*dd_sigma+(n**2+l**2-4d0*n**2*l**2)*dd_pi &
                      +(m**2+n**2*l**2)*dd_delta !zx,zx
    E2c_int(7,8,i)  = 1.5d0*n*l*(l**2-m**2)*dd_sigma &
                      +n*l*(1d0-2d0*(l**2-m**2))*dd_pi &
                      -n*l*(1d0-0.5d0*(l**2-m**2))*dd_delta !zx,x2-y2
    E2c_int(7,9,i)  = sqrt(3d0)*l*n*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                      +sqrt(3d0)*l*n*(l**2+m**2-n**2)*dd_pi &
                      -0.5d0*sqrt(3d0)*l*n*(l**2+m**2)*dd_delta !zx,3z2-r2
    E2c_int(7,10,i) = sqrt(3d0)*n+l*s2cda_sigma !zx,s2

! anion dx2-y2-orbital
    E2c_int(8,1,i)  = 0.5d0*sqrt(3d0)*(l**2-m**2)*scda_sigma !x2-y2,s
    E2c_int(8,2,i)  = (-1d0)*(0.5d0*sqrt(3d0)*l*(l**2-m**2)*pcda_sigma &
                      +l*(1d0-l**2+m**2)*pcda_pi) !x2-y2,x
    E2c_int(8,3,i)  = (-1d0)*(0.5d0*sqrt(3d0)*m*(l**2-m**2)*pcda_sigma &
                      -m*(1d0+l**2-m**2)*pcda_pi) !x2-y2,y
    E2c_int(8,4,i)  = (-1d0)*(0.5d0*sqrt(3d0)*n*(l**2-m**2)*pcda_sigma &
                      -n*(l**2-m**2)*pcda_pi) !x2-y2,z
    E2c_int(8,5,i)  = 1.5d0*l*m*(l**2-m**2)*dd_sigma+2d0*l*m*(m**2-l**2)*dd_pi &
                      0.5d0*l*m*(l**2-m**2)*dd_delta !x2-y2,xy
    E2c_int(8,6,i)  = 1.5d0*m*n*(l**2-m**2)*dd_sigma-m*n*(1d0+2d0*(l**2-m**2))*dd_pi &
                      +m*n*(1d0+0.5d0*(l**2-m**2))*dd_delta !x2-y2,yz
    E2c_int(8,7,i)  = 1.5d0*n*l*(l**2-m**2)*dd_sigma+n*l*(1d0-2d0*(l**2-m**2))*dd_pi &
                      -n*l*(1d0-0.5d0*(l**2-m**2))*dd_delta !x2-y2,zx
    E2c_int(8,8,i)  = 0.75d0*(l**2-m**2)**2*dd_sigma &
                      +(l**2+m**2-(l**2-m**2)**2)*dd_pi &
                      +(n**2+0.25d0*(l**2-m**2)**2)*dd_delta !x2-y2,x2-y2
    E2c_int(8,9,i)  = 0.5d0*sqrt(3d0)*(l**2-m**2)*(n**2-0.5d0*(l**2+m**2))*dd_sigma &
                      sqrt(3d0)*n**2*(m**2-l**2)*dd_pi &
                      +0.25d0*sqrt(3d0)*(1d0+n**2)*(l**2-m**2)*dd_delta !x2-y2,3z2-r2
    E2c_int(8,10,i) = 0.5d0*sqrt(3d0)*(l**2-m**2)*s2cda_sigma !x2-y2,s2

  end do

end subroutine calc_two_center_integral
