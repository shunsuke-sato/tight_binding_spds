program main
  use parallel
  use inputoutput
  use electronic_system
  implicit none

  call init_parallel
  call init_inputoutput

  call initialize_electronic_system
  call calc_bandstructure_zincblende

  call fin_inputoutput
  call fin_parallel

end program main
