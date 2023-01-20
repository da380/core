program test_SEM_1D

  use module_constants
  use module_LAPACK
  use module_SEM_1D
  implicit none

  integer(i4b) :: ngll,ndim,io,inode,ispec,i
  integer(i4b), dimension(:,:), allocatable :: ibool
  real(dp) :: x1,x2,dx,xs,u,sig,amp,start,finish
  type(mesh_1D) :: mesh
  class(rm), allocatable :: a
  type(rm) :: b

  ! build the mesh and boolean array
  ngll = 5
  x1 = 0.0_dp
  x2 = 1.0_dp
  dx = 0.001_dp
  mesh = build_mesh_1D(ngll,x1,x2,dx)
  call mesh%set_dirichlet()
  call build_boolean_scalar_1D(mesh,ibool,ndim)  

  ! build the system matrix
!  allocate(rm::a)
!  allocate(brm::a)
  allocate(psbrm::a)
  select type(a)
  class is(brm)
     call a%band(ngll-1,ngll-1)
  class is(psbrm)
     call a%band(ngll-1)
  end select
  call a%allocate(ndim,ndim)
  call build_laplace_matrix_1D(mesh,ibool,a)
  call cpu_time(start)
  call a%fac()
  call cpu_time(finish)
  print *, finish-start


  ! build the force and solve
  call b%allocate(ndim,1)
  xs = 0.3_dp
  sig = 0.05_dp
  amp = 1.0_dp
!  call build_delta_force_1D(mesh,ibool,xs,b)
  call build_gaussian_force_1D(mesh,ibool,xs,sig,amp,b)
  call a%bsub(b)
  

  ! write the solution
  open(newunit = io,file='test_SEM_1D.out')
  do ispec = 1,mesh%nspec
     do inode = 1,mesh%ngll
        i = ibool(inode,ispec)
        if(i /= 0) then
           u = b%elem(i,1)
        else
           u = 0.0_dp
        end if
        write(io,*) mesh%x(inode,ispec),u
     end do
  end do
  close(io)
  


end program test_SEM_1D




