subroutine calc_bandstructure_zincblende
  use global_variables
  implicit none
  integer :: ik, nkpoint_8
  real(8):: kvec_L(3), kvec_Gamma(3), kvec_X(3)
  real(8):: kvec_W(3), kvec_K(3)
  real(8),allocatable :: kk(:)
!LAPACK
  integer :: ndim
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

  ndim = 40
  lwork = 4*ndim**2+4*ndim+128
  allocate(work_lp(lwork),rwork(3*ndim-2),w(ndim))


!LAPACK
  allocate(kk(nkpoint)); kk = 0d0

  kvec_L=2d0*pi/lattice_const*0.5d0
  kvec_Gamma = 0d0
  kvec_X=0d0; kvec_X(1)=2d0*pi/lattice_const
  kvec_W=0d0; kvec_W(1)=2d0*pi/lattice_const; kvec_W(2)=2d0*pi/lattice_const*0.5d0
  kvec_K=2d0*pi/lattice_const*0.75d0; kvec_K(3)=0d0


  kvec0 = 0d0
  
  nkpoint_8 = nkpoint/8
  if(mod(nkpoint,8)/=0)then
    write(*,*)"mod(nkpoint,8)/=0"
    stop
  end if

!L to Gamma
  do ik = 1,nkpoint_8
    kvec(:,ik+0*nkpoint_8) = kvec_L(:) +(ik-1)*(kvec_Gamma(:)-kvec_L(:))/nkpoint_8
    if(ik /= 0)kk(ik+0*nkpoint_8)=kk(ik+0*nkpoint_8-1) &
        +sqrt(sum((kvec_Gamma(:)-kvec_L(:))**2))/nkpoint_8
  end do

!Gamma to X
  do ik = 1,nkpoint_8
    kvec(:,ik+1*nkpoint_8) = kvec_Gamma(:) +(ik-1)*(kvec_X(:)-kvec_Gamma(:))/nkpoint_8
    if(ik==0)then
      kk(ik+1*nkpoint_8)=kk(ik+1*nkpoint_8-1)+sqrt(sum((kvec_Gamma(:)-kvec_L(:))**2))/nkpoint_8
    else
      kk(ik+1*nkpoint_8)=kk(ik+1*nkpoint_8-1)+sqrt(sum((kvec_X(:)-kvec_Gamma(:))**2))/nkpoint_8
    end if
  end do

!X to W
  do ik = 1,nkpoint_8
    kvec(:,ik+2*nkpoint_8) = kvec_X(:) +(ik-1)*(kvec_W(:)-kvec_X(:))/nkpoint_8
    if(ik==0)then
      kk(ik+2*nkpoint_8)=kk(ik+2*nkpoint_8-1)+sqrt(sum((kvec_X(:)-kvec_Gamma(:))**2))/nkpoint_8
    else
      kk(ik+2*nkpoint_8)=kk(ik+2*nkpoint_8-1)+sqrt(sum((kvec_W(:)-kvec_X(:))**2))/nkpoint_8
    end if
  end do

!W to K
  do ik = 1,nkpoint_8
    kvec(:,ik+3*nkpoint_8) = kvec_W(:) +(ik-1)*(kvec_K(:)-kvec_W(:))/nkpoint_8
    if(ik==0)then
      kk(ik+3*nkpoint_8)=kk(ik+3*nkpoint_8-1)+sqrt(sum((kvec_W(:)-kvec_X(:))**2))/nkpoint_8
    else
      kk(ik+3*nkpoint_8)=kk(ik+3*nkpoint_8-1)+sqrt(sum((kvec_K(:)-kvec_W(:))**2))/nkpoint_8
    end if
  end do

!K to L
  do ik = 1,nkpoint_8
    kvec(:,ik+4*nkpoint_8) = kvec_K(:) +(ik-1)*(kvec_L(:)-kvec_K(:))/nkpoint_8
    if(ik==0)then
      kk(ik+4*nkpoint_8)=kk(ik+4*nkpoint_8-1)+sqrt(sum((kvec_K(:)-kvec_W(:))**2))/nkpoint_8
    else
      kk(ik+4*nkpoint_8)=kk(ik+4*nkpoint_8-1)+sqrt(sum((kvec_L(:)-kvec_K(:))**2))/nkpoint_8
    end if
  end do

!L to W
  do ik = 1,nkpoint_8
    kvec(:,ik+5*nkpoint_8) = kvec_L(:) +(ik-1)*(kvec_W(:)-kvec_L(:))/nkpoint_8
    if(ik==0)then
      kk(ik+5*nkpoint_8)=kk(ik+5*nkpoint_8-1)+sqrt(sum((kvec_L(:)-kvec_K(:))**2))/nkpoint_8
    else
      kk(ik+5*nkpoint_8)=kk(ik+5*nkpoint_8-1)+sqrt(sum((kvec_W(:)-kvec_L(:))**2))/nkpoint_8
    end if
  end do

!W to X
  do ik = 1,nkpoint_8
    kvec(:,ik+6*nkpoint_8) = kvec_W(:) +(ik-1)*(kvec_X(:)-kvec_W(:))/nkpoint_8
    if(ik==0)then
      kk(ik+6*nkpoint_8)=kk(ik+6*nkpoint_8-1)+sqrt(sum((kvec_W(:)-kvec_L(:))**2))/nkpoint_8
    else
      kk(ik+6*nkpoint_8)=kk(ik+6*nkpoint_8-1)+sqrt(sum((kvec_X(:)-kvec_W(:))**2))/nkpoint_8
    end if
  end do


!K to Gamma
  do ik = 1,nkpoint_8
    kvec(:,ik+7*nkpoint_8) = kvec_K(:) +(ik-1)*(kvec_Gamma(:)-kvec_K(:))/(nkpoint_8-1)
    if(ik==0)then
      kk(ik+7*nkpoint_8)=kk(ik+7*nkpoint_8-1)+sqrt(sum((kvec_X(:)-kvec_W(:))**2))/nkpoint_8
    else
      kk(ik+7*nkpoint_8)=kk(ik+7*nkpoint_8-1)+sqrt(sum((kvec_Gamma(:)-kvec_K(:))**2))/(nkpoint_8-1)
    end if
  end do


  call calc_two_center_integral
  call calc_zham_mat


  open(20,file="bandstructure.out")
  do ik = 1, nkpoint

    call zheev('V', 'U', ndim, zham_mat(:,:,ik), ndim, w, work_lp, lwork, rwork, info)
    write(20,"(999e26.16e3)")kk(ik)/kk(nkpoint),kvec(:,ik),w

  end do
  close(20)

end subroutine calc_bandstructure_zincblende
