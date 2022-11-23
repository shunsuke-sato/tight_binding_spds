subroutine input
  use global_variables
  implicit none

  nband = 2*2*(1+3+5+1)
  nkx = 9*2
  nky = 9*2
  nkz = 9*2
  nkpoint = nkx*nky*nkz
  

  lattice_const = 5.4635d0*angstrom

!  lattice_const = 100d0
  


end subroutine input
