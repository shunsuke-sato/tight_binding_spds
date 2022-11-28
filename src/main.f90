program main
  use parallel
  use inputoutput
  use electronic_system
  use electron_dynamics
  implicit none

  call init_parallel
  call init_inputoutput

  call initialize_electronic_system
!  call calc_bandstructure_zincblende

  call electron_dynamics_calculation

  call fin_inputoutput
  call fin_parallel

end program main
