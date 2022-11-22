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

end subroutine calc_zham_mat
