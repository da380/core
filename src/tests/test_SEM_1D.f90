program test_SEM_1D

  use module_constants
  use module_linalg
  use module_SEM_1D
  implicit none

  integer(i4b) :: k
  real(dp), dimension(3,3) :: ab
  type(matrix_list) :: al

  ab = 1.0_dp
  
  al = matrix_list(10,10)
  call al%print()
  call al%add(4,3,2,3,ab)
  call al%add(1,2,2.0_dp)
  call al%print()

end program test_SEM_1D




