program test_SEM_1D

  use module_constants
  use module_linalg
  use module_SEM_1D
  implicit none

  integer(i4b) :: ngll,n,kl,ku,ld,nn,info,io,inode,ispec,i
  integer(i4b), dimension(:), allocatable :: ipiv
  integer(i4b), dimension(:,:), allocatable :: ibool
  
  real(dp) :: x,x1,x2,dx,om
  real(dp), dimension(:,:), allocatable :: mm,kk,ss,ff

  type(mesh_1D) :: mesh
  type(matrix_list) :: list
  
  ! build the mesh
  ngll = 5
  x1 = 0.0_dp
  x2 = 1.0_dp
  dx = 0.001_dp
  mesh = mesh_1D(ngll,x1,x2,dx)
  call build_boolean_scalar_1D(mesh,ibool,n)

  ! allocate the banded matrices
  kl = ngll-1
  ku = ngll-1  
  ld = 2*kl+ku+1
  allocate(mm(ld,n),kk(ld,n),ss(ld,n))
  
  ! build the mass matrix
  list = matrix_list( n, n, herm = .true.)
  call build_laplace_mass_matrix_1D(mesh,ibool,list)
  call list%banded(kl,ku,mm)

  ! build stiffness matrix
  list = matrix_list( n, n, herm = .true.)
  call build_laplace_stiffness_matrix_1D(mesh,ibool,list)
  call list%banded(kl,ku,kk)  

  ! delete the matrix list
  call list%deallocate()

  ! form the combined matrix
  om = 10.0_dp
  ss = -om*om*mm + kk

  ! factor the matrix
  allocate(ipiv(n))
  call dgbtrf(n,n,kl,ku,ss,ld,ipiv,info)	
  call check(info == 0,'test_SEM_1D','problem with matrix factorisation')
  
  ! set the force term
  allocate(ff(n,1))
  ff = 0.0_dp
  ff(1,1) = 1.0_dp

  ! solve the linear system
  call dgbtrs('N',n,kl,ku,1,ss,ld,ipiv,ff,n,info) 
  call check(info == 0,'test_SEM_1D','problem with back-substitution')


  open(newunit = io,file = 'test_SEM_1D.out')
  do ispec = 1,mesh%nspec
     do inode = 1,mesh%ngll
        i = ibool(inode,ispec)
        write(io,*) mesh%x(inode,ispec),ff(i,1)
     end do
  end do
  close(io)
  
end program test_SEM_1D




