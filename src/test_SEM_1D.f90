program test_SEM_1D

  use module_constants
  use module_SEM_1D
  implicit none

  integer(i4b) :: ngll,nsec,inode,ispec
  real(dp) :: x1,x2,dx
  type(mesh_1D) :: mesh

  ngll = 5
  x1 = 0.0_dp
  x2 = 1.0_dp
  dx = 0.1_dp
  mesh = mesh_1D(ngll,x1,x2,dx)
  
  do ispec = 1,mesh%nspec
     do inode = 1,mesh%ngll
        print *, ispec,inode,mesh%x(inode,ispec)
     end do
  end do
  
end program test_SEM_1D



