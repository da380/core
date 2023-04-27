program test_interp

  
  use module_constants
  use module_interp
  use module_util
  implicit none

  integer(i4b), parameter :: n = 100, m = 700
  integer(i4b) :: i
  real(dp) :: x1,x2,dx,xx
  real(dp), dimension(n) :: x,f,fpp
  type(interp_1D_cubic) :: g

  
  x1 = 0.0_dp
  x2 = 1.0_dp
  dx = (x2-x1)/(n-1)
  xx = x1
  do i = 1,n
     x(i) = xx
     f(i) = sin(xx)
     xx = xx + dx
  end do

  xx = -0.001_dp
  i = 20
  i = hunt_list(x,xx,i)
  print *, i


  
  
end program test_interp
