subroutine calc_zham_mat
  use global_variables
  implicit  none
  integer :: ik,i,ia,ic
  complex(8) :: zphase


  zham_mat = 0d0
  do ik = 1, nkpoint

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
