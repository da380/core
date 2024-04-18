program test_random_field

  use module_constants
  use module_util
  use module_spherical_harmonics
  use module_random_fields
  implicit none

  integer(i4b) :: lmax,l,ith,iph,io,inode,ispec,i,nx,n
  real(dp) :: lambda,s,th,ph,x1,x2,x,y,sigma,dx,x0,th0,ph0
  real(dp), dimension(:), allocatable :: xx,ff
  real(dp), dimension(:,:), allocatable :: u
  complex(dpc), dimension(:), allocatable :: ulm
  type(gauss_legendre_grid) :: grid
  class(GRF_1D), allocatable :: fun

  class(GRF_S2), allocatable :: fun_S2
  

  !---------------------------------------------!
  !     test random functions on an interval    !
  !---------------------------------------------!
  
  x1 = 0.0_dp
  x2 = 1.0_dp
  if(.not.found_command_argument('-s',s)) stop 's missing'
  if(.not.found_command_argument('-lambda',lambda)) stop 'lambda missing'
  if(.not.found_command_argument('-sigma',sigma)) stop 'sigma missing'


  fun = GRF_1D_Fourier(x1,x2,lambda,s,sigma)
  
  call fun%sample(100)
  
  nx = 50*(x2-x1)/lambda
  dx = (x2-x1)/(nx-1)
  x0 = 0.5_dp*(x1+x2)
  open(newunit = io,file='random.out')
  do i = 1,nx
     x = x1 + (i-1)*dx
     write(io,*) x,fun%eval(1,x),fun%corr_eval(x,x0)
  end do
  close(io)
  
  
 
  
  
  

  
  

  
end program test_random_field
