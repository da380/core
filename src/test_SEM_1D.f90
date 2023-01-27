  program test_SEM_1D

  use module_constants
  use module_SEM_1D
  use module_LAPACK
  implicit none

  integer(i4b) :: ngll,n,kd,io,inode,ispec,i,l,j,k
  integer(i4b), dimension(:,:), allocatable :: ibool
  integer(i4b), dimension(:), allocatable :: ipiv
  real(dp) :: x,x1,x2,dx,xs,u,sig,amp,lambda
  real(dp), dimension(:), allocatable :: b
  real(dp), dimension(:,:), allocatable :: a
  type(mesh_1D) :: mesh



  ! build the mesh and boolean array
  ngll = 5
  x1 = 0.0_dp
  x2 = 1.0_dp
  dx = 0.001_dp
  mesh = build_mesh_1D(ngll,x1,x2,dx)
  call mesh%set_dirichlet()
  call build_boolean_scalar_1D(mesh,ibool,n=n)

  
  ! build the system matrix
  kd = ngll-1  
  call allocate_matrix_bhp(n,kd,a)
  call build_laplace_stiffness_matrix_1D_bhp(mesh,ibool,n,kd,a)

  ! factorise the matrix
  call factorise_matrix_bhp(n,kd,a)
  
  ! build the force and solve
  xs = 0.4_dp
  sig = 0.01_dp
  amp = 1.0_dp
  call allocate_vector(n,b)
  call build_gaussian_force_1D(mesh,ibool,xs,sig,amp,b)
  call backsub_matrix_bhp(n,kd,a,b)
  
 ! write the solution
  open(newunit = io,file='test_SEM_1D.out')
  do ispec = 1,mesh%nspec
     do inode = 1,mesh%ngll
        i = ibool(inode,ispec)
        if(i /= 0) then
           u = b(i)
        else
           u = 0.0_dp
        end if
        write(io,*) mesh%x(inode,ispec),u
     end do
       end do
  close(io)
  


end program test_SEM_1D




