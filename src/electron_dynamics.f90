module electron_dynamics
  use mpi
  use parallel
  use communication
  use math
  use constants
  use inputoutput
  use electronic_system
  use laser
  implicit none
  private

  public :: electron_dynamics_calculation, &
            initialize_electron_dynamics

  integer :: n_time_step
  real(8) :: dt, time_propagation
  real(8),allocatable :: tt(:)


contains
!----------------------------------------------------------------------------------------
  subroutine electron_dynamics_calculation
    implicit  none
    integer :: it
    real(8) :: jt_t(3), Act_t(3), num_elec

    call input_parameter_for_time_propagation
    call initialize_electron_dynamics



    call calc_vector_potential_time(tt(0), Act_t)
    call calc_current(Act_t, jt_t)
    call calc_num_electron(num_elec)

    if(if_root_global)then
      open(40,file="act_jt.out")
      write(40,"(A)")"# tt (a.u.), num_elec, Act(1:3) (a.u.), jt(1:3) (a.u.)"
      write(40,"(999e26.16e3)")tt(0),num_elec, Act_t, jt_t
    end if

    do it = 1, n_time_step
      if(if_root_global)write(*,*)'it=',it
      call dt_evolve_mod(it)

      call calc_vector_potential_time(tt(it), Act_t)
      call calc_current(Act_t, jt_t)
      call calc_num_electron(num_elec)
      if(if_root_global)then
        write(40,"(999e26.16e3)")tt(it),num_elec, Act_t, jt_t
      end if
    end do

    if(if_root_global)close(40)

  end subroutine electron_dynamics_calculation
!----------------------------------------------------------------------------------------
  subroutine initialize_electron_dynamics
    implicit none

    call set_equilibrium_density_matrix
    call init_laser    

  end subroutine initialize_electron_dynamics
!----------------------------------------------------------------------------------------
  subroutine input_parameter_for_time_propagation
    implicit none
    real(8) :: time_propagation_fs
    integer :: it

    call read_basic_input('n_time_step',n_time_step,val_default = -1)
    call read_basic_input('dt',dt,val_default = -1d0)
    call read_basic_input('time_propagation_fs',time_propagation_fs,val_default = -1d0)
    time_propagation = time_propagation_fs*fs

    if(n_time_step <= -1)then
      n_time_step = aint(time_propagation/dt)+1
      if(n_time_step <=0)n_time_step = 1
      if(if_root_global)write(*,"(A,2x,e26.16e3)")"Original dt (a.u.)=",dt

      dt = time_propagation/n_time_step
      if(if_root_global)write(*,"(A,2x,e26.16e3)")"Refined dt (a.u.)=",dt
    end if

    if(if_root_global)write(*,"(A,2x,I9)")"n_time_step =",n_time_step

    allocate(tt(0:n_time_step))
    do it = 0, n_time_step
      tt(it) = dt*it
    end do

  end subroutine input_parameter_for_time_propagation
!----------------------------------------------------------------------------------------
  subroutine  dt_evolve(it)
    implicit none
    integer,intent(in) :: it
    real(8) :: ttt, Act(3)
    Act = 0d0
    
!! === Start: propagation from tt(it-1) to tt(it-1)+dt/2 ===
    ttt = tt(it-1)
    call calc_vector_potential_time(ttt, Act)

    call dt_evolve_elec_system(Act,dt*0.5d0)
!! === End: propagation from tt(it-1) to tt(it-1)+dt/2 ===

!! === Start: propagation from tt(it-1)+dt/2 to tt(it) ===
    ttt = tt(it)
    call calc_vector_potential_time(ttt, Act)

    call dt_evolve_elec_system(Act,dt*0.5d0)
!! === End: propagation from tt(it-1)+dt/2 to tt(it) ===
    
  end subroutine dt_evolve
!----------------------------------------------------------------------------------------
  subroutine  dt_evolve_mod(it)
    implicit none
    integer,intent(in) :: it
    real(8) :: ttt, Act_1(3), Act_2(3)

    ttt = tt(it-1)
    call calc_vector_potential_time(ttt, Act_1)
    ttt = tt(it)
    call calc_vector_potential_time(ttt, Act_2)

    call dt_evolve_elec_system_mod(Act_1, Act_2 ,dt)

  end subroutine dt_evolve_mod
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
end module electron_dynamics
