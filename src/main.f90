program main
  use global_variables
  implicit none

  call init_parallel
  call init_input



  call fin_input
  call fin_parallel

end program main
