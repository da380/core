program test

  use module_constants
  use module_util

  integer(i4b), parameter :: n = 4,m = 100
  integer(i4b) :: i
  real(dp) :: x1,x2,dx,x,p,dpdx
  real(dp), dimension(n) :: coef

  coef(1) = 1.0_dp
  coef(2) = 1.0_dp
  coef(3) = 1.0_dp
  coef(4) = -3.0_dp
  
  x1 = 0.0_dp
  x2 = 1.0_dp
  dx = (x2-x1)/(m-1)

  x = x1
  do i = 1,m
     p  = poly(n,coef,x)
     dpdx = dpoly(n,coef,x)     
     print *, x,p - (1 + x + x*x - 3*x*x*x),dpdx - (1 + 2*x - 9*x*x)
     x = x+dx
  end do
  
  
end program test
