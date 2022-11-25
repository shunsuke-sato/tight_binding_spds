module electron_dynamics
  use mpi
  use parallel
  use communication
  use math
  use constants
  use inputoutput
  use electronic_system
  implicit none
  private

  public :: electron_dynamics_calculation, &
            initialize_electron_dynamics


contains
!----------------------------------------------------------------------------------------
  subroutine electron_dynamics_calculation
    implicit  none

    call initialize_electron_dynamics

  end subroutine electron_dynamics_calculation
!----------------------------------------------------------------------------------------
  subroutine initialize_electron_dynamics
    implicit none

    call set_equilibrium_density_matrix
    

  end subroutine initialize_electron_dynamics
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
end module electron_dynamics
