subroutine input
  use global_variables
  implicit none

  nband = 2*2*(1+3+5+1)
  nkx = 8
  nky = 8
  nkz = 8
  nkpoint = nkx*nky*nkz
  

  lattice_const = 5.4635d0*angstrom
  


end subroutine input
