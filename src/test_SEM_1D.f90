program test_SEM_1D

  use module_constants
  use module_LAPACK
  use module_SEM_1D
  implicit none


  integer(i4b) :: m,n
  type(real_matrix) :: b
  type(LU_real_matrix) :: a


  a = LU_real_matrix(3)
  a%elem(1,1) =  2.0_dp
  a%elem(2,1) = -1.0_dp
  a%elem(2,2) =  3.0_dp
  a%elem(3,3) =  1.0_dp
  call a%factorise()


  b = real_matrix(3,1)
  b%elem(1,1) = 1.0_dp
  b%elem(2,1) = 2.0_dp
  
  call a%backsub(b)
  print *, b%elem

  

end program test_SEM_1D




