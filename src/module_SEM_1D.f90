module module_SEM_1D

  use module_constants
  use module_error
  use module_quadrature
  use module_special_functions
  implicit none

  type mesh_1D
     logical :: allocated = .false.
     integer(i4b) :: nspec = 0
     integer(i4b) :: ngll  = 0
     real(dp), dimension(:,:), allocatable :: x
     real(dp), dimension(:), allocatable :: jac
     real(dp), dimension(:), allocatable :: w
     real(dp), dimension(:,:), allocatable :: hp
   contains
     procedure :: delete => delete_mesh_1D
  end type mesh_1D

  interface mesh_1D
     procedure :: build_mesh_1D
     procedure :: build_mesh_1D_constant_dx
     procedure :: build_mesh_1D_simple
  end interface mesh_1D
  
contains

  subroutine delete_mesh_1D(self)
    class(mesh_1D), intent(inout) :: self
    if(.not.self%allocated) return
    self%nspec = 0
    self%ngll  = 0
    deallocate(self%x)    
    deallocate(self%jac)
    deallocate(self%w)
    deallocate(self%hp)
    self%allocated = .false.
    return
  end subroutine delete_mesh_1D

  type(mesh_1D) function build_mesh_1D(ngll,x,dx) result(mesh)
    integer(i4b), intent(in) :: ngll
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in) :: dx

    integer(i4b) :: ispec,nspec,isec,nspec_loc,jspec,inode,nsec
    real(dp) :: x1,x2,dx_loc,x1_loc,x2_loc
    real(dp), dimension(ngll) :: h
    type(gauss_lobatto_quadrature) :: quad


    nsec = size(dx)
    call check(nsec+1 == size(x),'build_mesh_1D','input arrays have wrong dimensions')
    
    ! work out number of spectral elements
    nspec = 0
    do isec = 1,nsec
       x1 = x(isec)
       x2 = x(isec+1)
       nspec_loc = (x2-x1)/dx(isec)
       nspec_loc = max(nspec_loc,1)
       nspec = nspec + nspec_loc
    end do


    ! set up the mesh arrays
    mesh%nspec = nspec
    mesh%ngll = ngll
    allocate(mesh%x(ngll,nspec))
    allocate(mesh%jac(nspec))
    allocate(mesh%hp(ngll,ngll))
    mesh%allocated = .true.

    ! set the quadrature scheme
    call quad%set(ngll)
    mesh%w = quad%w

    ! store the lagrange derivatives
    do inode = 1,ngll
       call lagrange_polynomial(quad%x(inode),ngll,quad%x,h,mesh%hp(inode,:))
    end do

    
    ! build up the mesh
    ispec = 0
    do isec = 1,nsec
       x1 = x(isec)
       x2 = x(isec+1)
       nspec_loc = (x2-x1)/dx(isec)
       nspec_loc = max(nspec_loc,1)
       dx_loc = (x2-x1)/(nspec_loc)
       do jspec = 1,nspec_loc
          ispec = ispec+1
          x1_loc = x1+(jspec-1)*dx_loc
          x2_loc = x1_loc+dx_loc
          do inode = 1,ngll
             mesh%x(inode,ispec) = x1_loc + 0.5_dp*dx_loc*(quad%x(inode)+1.0_dp)
          end do
       end do
    end do
        
    return
  end function build_mesh_1D
  
  type(mesh_1D) function build_mesh_1D_constant_dx(ngll,x,dx) result(mesh)
    integer(i4b), intent(in) :: ngll
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: dx
    real(dp), dimension(:), allocatable :: dx_loc
    integer(i4b) :: nsec
    nsec = size(x)-1
    allocate(dx_loc(nsec))
    dx_loc = dx
    mesh = build_mesh_1D(ngll,x,dx_loc)
    return
  end function build_mesh_1D_constant_dx

  type(mesh_1D) function build_mesh_1D_simple(ngll,x1,x2,dx) result(mesh)
    integer(i4b), intent(in) :: ngll
    real(dp), intent(in) :: x1,x2,dx
    real(dp), dimension(2) :: x_loc
    real(dp), dimension(1) :: dx_loc
    x_loc(1) = x1
    x_loc(2) = x2
    dx_loc(1) = dx
    mesh = build_mesh_1D(ngll,x_loc,dx_loc)
    return
  end function build_mesh_1D_simple


  
end module module_SEM_1D
