program main
  use parallel
  use inputoutput
  use electronic_system
  use electron_dynamics
#ifdef profile
  use profile_m
#endif
  implicit none

  call init_parallel
#ifdef profile
  call init_profile
#endif
  
  call init_inputoutput

  call initialize_electronic_system
!  call calc_bandstructure_zincblende

  call electron_dynamics_calculation

  call fin_inputoutput

#ifdef profile
  call fin_profile
#endif  
  call fin_parallel

end program main
