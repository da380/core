program test_random_field

  use module_constants
  use module_util
  use module_spherical_harmonics
  use module_random_fields
  implicit none

  integer(i4b) :: lmax,l,ith,iph,io,inode,ispec,i,nx,n
  real(dp) :: lambda,s,th,ph,x1,x2,x,y,sigma,dx,x0
  real(dp), dimension(:), allocatable :: xx,ff
  real(dp), dimension(:,:), allocatable :: v
  type(gauss_legendre_grid) :: grid
!  type(gaussain_random_scalar_field_sphere) :: u
  type(interp_1D_cubic) :: fun1,cfun1
  type(interp_1D_cubic) :: fun2,cfun2
  class(GRF_1D), allocatable :: rfun1,rfun2


  !---------------------------------------------!
  !     test random functions on an interval    !
  !---------------------------------------------!
  
  x1 = 0.0_dp
  x2 = 1.0_dp
  if(.not.found_command_argument('-s',s)) stop 's missing'
  if(.not.found_command_argument('-lambda',lambda)) stop 'lambda missing'
  if(.not.found_command_argument('-sigma',sigma)) stop 'sigma missing'


  rfun1 =     GRF_1D_SEM(x1,x2,lambda,s,sigma)
  rfun2 = GRF_1D_Fourier(x1,x2,lambda,s,sigma)

  x0 = 0.2_dp
  call rfun1%corr(x0,cfun1)
  call rfun1%realise(fun1)

  call rfun2%corr(x0,cfun2)
  call rfun2%realise(fun2)
  
  nx = 20*(x2-x1)/lambda
  dx = (x2-x1)/(nx-1)
  x0 = 0.3_dp
  open(newunit = io,file='random.out')
  do i = 1,nx
     x = x1 + (i-1)*dx
     write(io,*) x,cfun1%f(x),cfun2%f(x),cfun1%f(x)-cfun2%f(x)
  end do
  close(io)
  


  
 
  
  
  !----------------------------------------!
  !   test random functions on a sphere    !
  !----------------------------------------!

  ! get arguments
!  if(.not.found_command_argument('-lmax',lmax)) stop 'lmax missing'
!  if(.not.found_command_argument('-s',s)) stop 's missing'
!  if(.not.found_command_argument('-lambda',lambda)) stop 'lambda missing'
  
  ! build the grid
!  call grid%build(lmax)  

  ! allocate the random field
!  call u%build(grid,lambda,s)

  ! allocate coefficient array and spatial array
!  allocate(v(grid%nph,grid%nth))
  
  ! form realisation of the random field
!  call u%realise()

!  call grid%SH_itrans(u%ulm,v)
  
  ! write out the field
!  open(newunit = io,file='random.out')
!  write(io,*) grid%nth,grid%nph,0.0_dp
!  do ith = 1,grid%nth
!     th = grid%th(ith)
!     do iph = 1,grid%nph
!        ph = grid%ph(iph)
!        write(io,*) ph,th,v(iph,ith)
!     end do
!  end do
!  close(io)  
  
  
end program test_random_field
