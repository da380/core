program test_SEM_1D

  use module_constants
  use module_LAPACK
  use module_SEM_1D
  implicit none


  integer(i4b) :: m,n,i,j
  real(dp) :: start,end
  type(rmat) :: a,b,c

  a = rmat(3)
  b = rmat(3,1)
  do i = 1,a%m
     call a%inc(i,i,1.0_dp)
     call b%inc(i,1,1.0_dp)     
     do j = 1,a%n
        call a%inc(i,j,0.3_dp*(i+j))
     end do
  end do
  c = rmat(a)
  call a%LU()
  call a%LU_backsub(b)
  print *, matmul(c%elem,b%elem)


  

  


end program test_SEM_1D




