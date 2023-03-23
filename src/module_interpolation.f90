module module_interpolation

  use module_constants
  use module_error
  use module_util
  implicit none

  type, abstract :: interpolant_1D
     private
     logical :: lhunt = .true.     
     integer(i4b) :: it = 0
     real(dp), dimension(:), pointer :: x
     real(dp), dimension(:), pointer :: f
   contains
     procedure, public :: deallocate => deallocate_interpolant_1D
     procedure :: find => find_interpolant_1D
     procedure(evaluate_interpolant_1D), deferred, public:: evaluate
     procedure(derivative_interpolant_1D), deferred, public:: derivative
  end type interpolant_1D


  type, extends(interpolant_1D) :: linear_interpolant_1D
   contains
     procedure :: evaluate => evaluate_linear_interpolant_1D
     procedure :: derivative => derivative_linear_interpolant_1D
  end type linear_interpolant_1D

  type, extends(interpolant_1D) :: cubic_interpolant_1D
     real(dp), dimension(:), allocatable :: fp2
   contains
     procedure :: deallocate => deallocate_cubic_interpolant_1D
     procedure :: evaluate => evaluate_cubic_interpolant_1D
     procedure :: derivative => derivative_cubic_interpolant_1D
  end type cubic_interpolant_1D
  
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
     procedure :: set_linear_interpolant_1D
  end interface linear_interpolant_1D


  interface cubic_interpolant_1D
     procedure :: set_cubic_interpolant_1D
  end interface cubic_interpolant_1D
  
contains


  !==========================================================!
  !                 procedures for basic type                !
  !==========================================================!
  
  subroutine deallocate_interpolant_1D(self)
    class(interpolant_1D), intent(inout) :: self
    nullify(self%x)
    nullify(self%f)
    self%lhunt = .true.
    return
  end subroutine deallocate_interpolant_1D


  integer(i4b) function find_interpolant_1D(self,x) result(i)
    class(interpolant_1D), intent(inout) :: self
    real(dp), intent(in) :: x
    i = hunt_list(self%x,x,self%it)
    if(self%lhunt) self%it = i
    return
  end function find_interpolant_1D


  !==========================================================!
  !           procedures for linear interpolant_1Ds          !
  !==========================================================!
  
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


  real(dp) function derivative_linear_interpolant_1D(self,x) result(fp)
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
    fp = (f2-f1)/(x2-x1)
    return
  end function derivative_linear_interpolant_1D
  
  
  type(linear_interpolant_1D) function set_linear_interpolant_1D(x,f,hunt) result(g)
    real(dp), dimension(:), target, intent(in) :: x
    real(dp), dimension(:), target, intent(in) :: f
    logical, optional, intent(in) :: hunt
    call check(size(x) == size(f), 'set_linear_interpolant_1D','arrays of different lengths')
    g%x => x
    g%f => f
    if(present(hunt)) g%lhunt = hunt
    return
  end function set_linear_interpolant_1D


  !==========================================================!
  !        procedures for cubic spline interpolant_1Ds       !
  !==========================================================!

  subroutine deallocate_cubic_interpolant_1D(self)
    class(cubic_interpolant_1D), intent(inout) :: self
    deallocate(self%fp2)
    nullify(self%x)
    nullify(self%f)
    self%lhunt = .true.    
    return
  end subroutine deallocate_cubic_interpolant_1D
  

  real(dp) function evaluate_cubic_interpolant_1D(self,x) result(f)
    class(cubic_interpolant_1D), intent(inout) :: self
    real(dp), intent(in) :: x
    integer(i4b) :: i1,i2
    real(dp) :: a,b,h,x1,x2,f1,f2,f21,f22    
    i1 = self%find(x)
    i2 = i1+1
    x1  = self%x(i1)
    x2  = self%x(i2)
    f1  = self%f(i1)
    f2  = self%f(i2)
    f21 = self%fp2(i1)
    f22 = self%fp2(i2)
    h   = x2-x1
    a   = (x2-x)/h
    b   = (x-x1)/h
    f = a*f1 + b*f2 + ((a*a*a-a)*f21 +  (b*b*b-b)*f22)*(h*h)/6.0_dp
    return
  end function evaluate_cubic_interpolant_1D


  real(dp) function derivative_cubic_interpolant_1D(self,x) result(fp)
    class(cubic_interpolant_1D), intent(inout) :: self
    real(dp), intent(in) :: x
        integer(i4b) :: i1,i2
    real(dp) :: a,b,h,x1,x2,f1,f2,f21,f22
    i1 = self%find(x)
    i2 = i1+1
    x1  = self%x(i1)
    x2  = self%x(i2)
    f1  = self%f(i1)
    f2  = self%f(i2)
    f21 = self%fp2(i1)
    f22 = self%fp2(i2)
    h   = x2-x1
    a   = (x2-x)/h
    b   = (x-x1)/h
    fp = (f2-f1)/h  - (3.0_dp*a*a-1)*h*f21/6.0_dp  &
                    + (3.0_dp*b*b-1)*h*f22/6.0_dp
    return
  end function derivative_cubic_interpolant_1D

  
  type(cubic_interpolant_1D) function set_cubic_interpolant_1D(x,f,fpl,fpr,hunt) result(g)
    real(dp), dimension(:), target, intent(in) :: x
    real(dp), dimension(:), target, intent(in) :: f
    real(dp), optional, intent(in) :: fpl,fpr
    logical, optional, intent(in) :: hunt
    logical :: lnat,rnat
    integer(i4b) :: i,n,info
    real(dp), dimension(:), allocatable :: dl,d,du
    real(dp), dimension(:,:), allocatable :: b

    ! set basic data
    n = size(x)
    call check(n == size(f), 'set_linear_interpolant_1D','arrays of different lengths')
    g%x => x
    g%f => f
    allocate(g%fp2(n))
    g%fp2 = 0.0_dp

    ! deal with optional arguments
    if(present(fpl)) then
       lnat = .false.
    else
       lnat = .true.
    end if
    if(present(fpr)) then
       rnat = .false.
    else
       rnat = .true.
    end if
    if(present(hunt)) g%lhunt = hunt

    ! build the matrix
    allocate(dl(n-1),d(n),du(n-1),b(n,1))

    ! upper diagonal
    if(lnat) then
       du(1) = 0.0_dp
    else
       du(1) = (x(2)-x(1))/6.0_dp
    end if
    do i = 2,n-1
       du(i) = (x(i+1)-x(i))/6.0_dp
    end do
    
    ! lower diagonl
    do i = 2,n-1
       dl(i-1) = (x(i)-x(i-1))/6.0_dp       
    end do
    if(rnat) then
       dl(n-1) = 0.0_dp
    else
       dl(n-1) = (x(n)-x(n-1))/6.0_dp
    end if

    ! diagonal and right hand side
    do i = 2,n-1
       d(i)   =  (x(i+1)-x(i-1))/3.0_dp
       b(i,1) =  (f(i+1)-f(i))/(x(i+1)-x(i)) &
                -(f(i)-f(i-1))/(x(i)-x(i-1))
    end do
    if(lnat) then
       d(1) = 1.0_dp
       b(1,1) = 0.0_dp
    else
       d(1) = (x(2)-x(1))/3.0_dp
       b(1,1) = (f(2)-f(1))/(x(2)-x(1)) - fpl
    end if
    if(rnat) then
       d(n) = 1.0_dp
       b(n,1) = 0.0_dp
    else
       d(n) = -(x(n)-x(n-1))/3.0_dp
       b(1,1) = (f(n)-f(n-1))/(x(n)-x(n-1)) - fpr
    end if
    
    ! solve the tridiagonal system
    call dgtsv (n, 1, dl, d, du, b, n, info)
    g%fp2 = b(:,1)
    
    return
  end function set_cubic_interpolant_1D

  
end module module_interpolation
