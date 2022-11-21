program main
  use global_variables
  implicit none

  call init_parallel
  call init_input

  call input
  call initialize

  call calc_bandstructure_zincblende

  call fin_input
  call fin_parallel

end program main
