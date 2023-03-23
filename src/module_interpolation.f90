module module_interpolation

  use module_constants
  use module_util
  implicit none

  type, abstract :: interpolant_1D
     logical, public :: hunt = .true.
     integer(i4b), private :: it = 0
     real(dp), dimension(:), pointer, private :: x
     real(dp), dimension(:), pointer, private :: f
   contains
     procedure :: deallocate => deallocate_interpolant_1D
     procedure, private :: find => find_interpolant_1D
     procedure(evaluate_interpolant_1D), deferred:: evaluate
     procedure(derivative_interpolant_1D), deferred:: derivative
  end type interpolant_1D


  type, extends(interpolant_1D) :: linear_interpolant_1D
   contains
     procedure :: evaluate => evaluate_linear_interpolant_1D
     procedure :: derivative => derivative_linear_interpolant_1D
  end type linear_interpolant_1D
  
  abstract interface

     real(dp) function evaluate_interpolant_1D(self,x) result(f)
       use module_constants
       import interpolant_1D
       class(interpolant_1D), intent(inout) :: self
       real(dp), intent(in) :: x
     end function evaluate_interpolant_1D

     real(dp) function derivative_interpolant_1D(self,x) result(f)
       use module_constants
       import interpolant_1D
       class(interpolant_1D), intent(inout) :: self
       real(dp), intent(in) :: x
     end function derivative_interpolant_1D    
     
  end interface

  interface linear_interpolant_1D
     procedure :: build_linear_interpolant_1D
  end interface linear_interpolant_1D
  
contains

  subroutine deallocate_interpolant_1D(self)
    class(interpolant_1D), intent(inout) :: self
    nullify(self%x)
    nullify(self%f)
    return
  end subroutine deallocate_interpolant_1D

  integer(i4b) function find_interpolant_1D(self,x) result(i)
    class(interpolant_1D), intent(inout) :: self
    real(dp), intent(in) :: x
    i = hunt_list(self%x,x,self%it)
    if(self%hunt) self%it = i
    return
  end function find_interpolant_1D

  real(dp) function evaluate_linear_interpolant_1D(self,x) result(f)
    class(linear_interpolant_1D), intent(inout) :: self
    real(dp), intent(in) :: x
    integer(i4b) :: i1,i2
    real(dp) :: x1,x2,f1,f2,h
    i1 = self%find(x)
    i2 = i1+1
    x1 = self%x(i1)
    x2 = self%x(i2)    
    f1 = self%f(i1)
    f2 = self%f(i2)
    h = (x-x1)/(x2-x1)
    f = f1 + h*(f2-f1)
    return
  end function evaluate_linear_interpolant_1D


  real(dp) function derivative_linear_interpolant_1D(self,x) result(df)
    class(linear_interpolant_1D), intent(inout) :: self
    real(dp), intent(in) :: x
    integer(i4b) :: i1,i2
    real(dp) :: x1,x2,f1,f2,h
    i1 = self%find(x)
    i2 = i1+1
    x1 = self%x(i1)
    x2 = self%x(i2)    
    f1 = self%f(i1)
    f2 = self%f(i2)
    df = (f2-f1)/(x2-x1)
    return
  end function derivative_linear_interpolant_1D
  
  
  type(linear_interpolant_1D) function build_linear_interpolant_1D(x,f) result(g)
    real(dp), dimension(:), target, intent(in) :: x
    real(dp), dimension(:), target, intent(in) :: f
    g%x => x
    g%f => f
    return
  end function build_linear_interpolant_1D
  
  
end module module_interpolation
