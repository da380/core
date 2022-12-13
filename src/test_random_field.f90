program test_random_field

  use module_constants
  use module_util
  use module_spherical_harmonics
  use module_random_fields
  implicit none

  integer(i4b) :: lmax,l,ith,iph,io,inode,ispec,i,nx,n
  real(dp) :: lambda,s,th,ph,x1,x2,x,y,sigma,dx,x0,th0,ph0
  real(dp), dimension(:), allocatable :: xx,ff
  real(dp), dimension(:,:), allocatable :: v
  complex(dpc), dimension(:), allocatable :: vlm
  type(gauss_legendre_grid) :: grid
  type(interp_1D_cubic) :: fun,cfun
  class(GRF_1D), allocatable :: rfun

  class(GRF_S2), allocatable :: rfun_S2
  

  !---------------------------------------------!
  !     test random functions on an interval    !
  !---------------------------------------------!
  
  x1 = 0.0_dp
  x2 = 1.0_dp
  if(.not.found_command_argument('-s',s)) stop 's missing'
  if(.not.found_command_argument('-lambda',lambda)) stop 'lambda missing'
  if(.not.found_command_argument('-sigma',sigma)) stop 'sigma missing'


  rfun = GRF_1D_Fourier(x1,x2,lambda,s,sigma)

  x0 = 0.1_dp

  call rfun%corr(x0,cfun)
  call rfun%realise(fun)
  
  nx = 50*(x2-x1)/lambda
  dx = (x2-x1)/(nx-1)
  x0 = 0.3_dp
  open(newunit = io,file='random.out')
  do i = 1,nx
     x = x1 + (i-1)*dx
         write(io,*) x,fun%f(x),cfun%f(x)
  end do
  close(io)
  
  
 
  
  
  !----------------------------------------!
  !   test random functions on a sphere    !
  !----------------------------------------!


  rfun_S2 = build_GRF_S2_SH(lambda,s,sigma)
  
  ! build the grid
  lmax = max(128,rfun_S2%lmax)

  call grid%build(lmax)  

  ! allocate coefficient array and spatial array
  allocate(vlm(grid%ncoef_r))
  allocate(v(grid%nph,grid%nth))
  
  ! form realisation of the random field
  call rfun_S2%realise(vlm)
  th0 = -20.0_dp
  ph0 = 100.0_dp
  th0 = (90.0_dp-th0)*deg2rad
  ph0 = ph0*deg2rad
  call rfun_S2%corr(th0,ph0,vlm)
  

  ! transform to spatial field
  call grid%SH_itrans(vlm,v)
  
  ! write out the field
  open(newunit = io,file='random_S2.out')
  write(io,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        write(io,*) ph,th,v(iph,ith)
     end do
  end do
  close(io)  
  
  
end program test_random_field
