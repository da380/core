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
  
  x1 = 0.000_dp
  x2 = 1.0_dp
  if(.not.found_command_argument('-s',s)) stop 's missing'
  if(.not.found_command_argument('-lambda',lambda)) stop 'lambda missing'
  if(.not.found_command_argument('-sigma',sigma)) stop 'sigma missing'


  fun = GRF_1D_SEM(x1,x2,0,lambda,s,sigma)
!  fun = GRF_1D_Fourier(x1,x2,lambda,s,sigma)
  
  call fun%realise()
  
  nx = 50*(x2-x1)/lambda
  dx = (x2-x1)/(nx-1)
  x0 = 0.5_dp*(x1+x2)
  open(newunit = io,file='random.out')
  do i = 1,nx
     x = x1 + (i-1)*dx
     write(io,*) x,fun%eval(x),fun%corr_eval(x,x0)
  end do
  close(io)
  
  
 
  
  
  !----------------------------------------!
  !   test random functions on a sphere    !
  !----------------------------------------!

  ! generate the random field
!  fun_S2 = build_GRF_S2_SH(lambda,s,sigma,asph = .true.)
!  call fun_S2%realise()

  ! build a GL-grid
!  lmax = max(128,fun_S2%degree())
!  call grid%build(lmax,0)
!  allocate(ulm(grid%ncoef_r))
!  allocate(u(grid%nph,grid%nth))
  
!  th0 = -20.0_dp
!  ph0 = 100.0_dp
!  th0 = (90.0_dp-th0)*deg2rad
!  ph0 = ph0*deg2rad
!  call fun_S2%coef(lmax,ulm)
!  !  call fun_S2%corr_coef(th0,ph0,lmax,ulm)  
!  call grid%SH_itrans(ulm,u)
  
 
  ! write out the field
!  open(newunit = io,file='random_S2.out')
!  write(io,*) grid%nth,grid%nph,0.0_dp
!  do ith = 1,grid%nth
!     th = grid%th(ith)
!     do iph = 1,grid%nph
!        ph = grid%ph(iph)
!        write(io,*) ph,th,u(iph,ith)
!     end do
!  end do
!  close(io)    

  
  

  
end program test_random_field
