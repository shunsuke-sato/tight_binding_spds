subroutine input
  use global_variables
  implicit none

  nband = 2*2*(1+3+5+1)
  nkx = 4
  nky = 4
  nkz = 4
  nkpoint = nkx*nky*nkz
  

  lattice_const = 5.4635d0*angstrom
  


end subroutine input
