program test_SEM_1D

  use module_constants
  use module_linalg
  use module_SEM_1D
  implicit none

  integer(i4b), parameter :: m = 3,n = 3
  integer(i4b) :: k,i,j
  real(dp), dimension(m,n) :: a
  type(matrix_list) :: al


  
  al = matrix_list(m,n,herm = .true.)
  call al%add(1,1,1.0_dp)
  call al%add(1,2,0.5_dp)

  call al%dense(a)
  do i = 1,m
     print *, a(i,:)
  end do


end program test_SEM_1D




