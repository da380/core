module module_interp

  use module_constants
  use module_error
  use module_util
  implicit none



  type :: interp_1D
     private
     logical, public :: allocated = .false.
     integer(i4b) :: n
     integer(i4b) :: i_save = 0
     real(dp), dimension(:), allocatable :: xx
     real(dp), dimension(:), allocatable :: ff
     real(dp) :: x_save,dx_save
     real(dp) :: fp_save
   contains
     procedure :: set => set_interp_1D
     procedure :: deallocate => deallocate_interp_1D
     procedure :: find   => find_interp_1D
     procedure :: f    => value_interp_1D
     procedure :: fp => derivative_interp_1D
  end type interp_1D
  private :: set_interp_1D,find_interp_1D, &
             value_interp_1D,deallocate_interp_1D
  
  
  type, extends(interp_1D) :: interp_1D_cubic
     private
     logical :: lnat = .true.
     logical :: rnat = .true.
     real(dp) :: fpl,fpr
     real(dp), dimension(:), allocatable :: ff2
     real(dp) :: fpp_save
   contains
     procedure :: parms => parameters_cubic_spline
     procedure :: set => set_interp_1D_cubic
     procedure :: deallocate => deallocate_interp_1D_cubic
     procedure :: f => value_interp_1D_cubic
     procedure :: fp => derivative_interp_1D_cubic
     procedure :: fpp => derivative2_interp_1D_cubic
  end type interp_1D_cubic
 private :: set_interp_1D_cubic,value_interp_1D_cubic, &
             deallocate_interp_1D_cubic

 
 type :: interp_2D
    private
    logical :: allocated = .false.
    integer(i4b) :: nx
    integer(i4b) :: ny
    real(dp), dimension(:), allocatable :: xx
    real(dp), dimension(:), allocatable :: yy
    real(dp), dimension(:,:), allocatable :: ff
    integer(i4b) :: ix_save = 0
    integer(i4b) :: iy_save = 0
    real(dp) :: x_save,dx_save
    real(dp) :: y_save,dy_save
  contains
    procedure :: set => set_interp_2D
    procedure :: deallocate => deallocate_interp_2D
    procedure :: find   => find_interp_2D
    procedure :: f    => value_interp_2D    
 end type interp_2D


  type, extends(interp_2D) :: interp_2D_bicubic_spline
     type(interp_1D_cubic), dimension(:), allocatable :: cs
   contains
     procedure :: set => set_interp_2D_bicubic_spline
     procedure :: deallocate => deallocate_interp_2D_bicubic_spline
     procedure :: f    => value_interp_2D_bicubic_spline    
 end type interp_2D_bicubic_spline

 
contains


  !======================================================================!
  !             general procedures for interpolation in 1D               !
  !======================================================================!
  
  
  subroutine set_interp_1D(self,xx,ff)
    implicit none
    class(interp_1D) :: self
    real(dp), dimension(:), intent(in) :: xx,ff
    
    integer(i4b) :: n,j,i
    
    ! check the size of the arrays
    n = size(xx)
    call check(n == size(ff),'set_linear_interp_1D', &
                             'dimensions of arrays do not match')
    
    ! check the arrays are consecutive
    j = 0
    do i = 1,n-1
       if(xx(i+1) >  xx(i)) then
          j = j+1
       else if(xx(i+1) < xx(i)) then
          j = j-1
       end if
    end do
    call check(j == n-1 .or. j == -n+1,'set_linear_interp_1D','array not consecutive')

    
    ! store the dimensions
    self%n = n

    ! store the data
    self%xx = xx
    self%ff = ff

    ! set dx_save
    self%dx_save = abs(xx(n)-xx(1))/5.0_dp

    ! build done
    self%allocated = .true.
    
    return
  end subroutine set_interp_1D

  subroutine find_interp_1D(self,x,i)
    use module_util, only: bisect_list,hunt_list
    implicit none
    class(interp_1D) :: self
    real(dp), intent(in) :: x
    integer(i4b), intent(out) :: i

    integer(i4b) :: il,iu,im,n

    ! is it worth hunting? 
    if(self%i_save /= 0 .and. abs(x-self%x_save) > self%dx_save) then
       self%i_save = 0
    end if

    ! find the point
    i = hunt_list(self%xx,x,self%i_save)

    ! _save the previous values
    self%i_save = i
    self%x_save = x
    
    return
  end subroutine find_interp_1D

  function value_interp_1D(self,x) result(f)
    implicit none
    class(interp_1D) :: self
    real(dp), intent(in) :: x

    integer(i4b) :: i
    real(dp) :: x1,x2,f1,f2,f
    
    call self%find(x,i)
    x1 = self%xx(i)
    x2 = self%xx(i+1)

    f1 = self%ff(i)
    f2 = self%ff(i+1)

    self%fp_save = (f2-f1)/(x2-x1)
    f = f1 + self%fp_save*(x-x1)
    
    return
  end function value_interp_1D


  function derivative_interp_1D(self,x) result(fp)
    implicit none
    class(interp_1D) :: self
    real(dp), intent(in) :: x

    integer(i4b) :: i
    real(dp) :: f,fp
    
    if(x /= self%x_save) f = self%f(x)
    fp = self%fp_save
    
    return
  end function derivative_interp_1D


  subroutine deallocate_interp_1D(self)
    implicit none
    class(interp_1D) :: self
    deallocate(self%xx,self%ff)
    self%allocated = .false.    
    return
  end subroutine deallocate_interp_1D
  
    
  !==================================================================!
  !              routines for cubic spline interpolation             !
  !==================================================================!

  subroutine build_cubic_spline(n,x,f,fpp,fpl,fpr)
    ! builds the cubic spline coefficients, fpp, given
    ! x and f values for a function. Derivative values
    ! at either endpoint can be optionally set, with
    ! the default being a natural spline    
    integer(i4b), intent(in) :: n
    real(dp), dimension(n), intent(in) :: x
    real(dp), dimension(n), intent(in) :: f
    real(dp), dimension(n), intent(out), target :: fpp
    real(dp), intent(in), optional :: fpl
    real(dp), intent(in), optional :: fpr

    logical :: lnat,rnat
    integer(i4b) :: i,info
    real(dp), dimension(n-1) :: dl,du
    real(dp), dimension(n) :: d
    real(dp), dimension(:,:), pointer :: b
    
    ! deal with optional arguments
    lnat = .not.present(fpl)
    rnat = .not.present(fpr)

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
       fpp(i) =  (f(i+1)-f(i))/(x(i+1)-x(i)) &
                -(f(i)-f(i-1))/(x(i)-x(i-1))
    end do    
    if(lnat) then
       d(1) = 1.0_dp
       fpp(1) = 0.0_dp
    else
       d(1) = (x(2)-x(1))/3.0_dp      
       fpp(1) = (f(2)-f(1))/(x(2)-x(1)) - fpl
    end if
    if(rnat) then
       d(n) = 1.0_dp
       fpp(n) = 0.0_dp
    else
       d(n) = -(x(n)-x(n-1))/3.0_dp
       fpp(n) = (f(n)-f(n-1))/(x(n)-x(n-1)) - fpr
    end if    
    
    ! solve the tridiagonal system
    b(1:n,1:1) => fpp
    call dgtsv (n, 1, dl, d, du, b, n, info)
    call check(info == 0,'build cubic spline','problem solving linear system')
    
    return
  end subroutine build_cubic_spline


  subroutine evaluate_cubic_spline(n,x,f,fpp,xx,i,ff,ffp,ffpp)
    ! evaluates the value of a cubic spline at a point xx
    ! can optionally return the value of the function or
    ! its first two derivatives. 
    integer(i4b), intent(in) :: n
    real(dp), dimension(n), intent(in) :: x
    real(dp), dimension(n), intent(in) :: f
    real(dp), dimension(n), intent(in) :: fpp
    real(dp), intent(in) :: xx
    integer(i4b), intent(inout), optional :: i
    real(dp), intent(out), optional :: ff
    real(dp), intent(out), optional :: ffp    
    real(dp), intent(out), optional :: ffpp

    integer(i4b) :: i1,i2
    real(dp) :: x1,x2,f1,f2,fpp1,fpp2,h,a,b

    ! bracket the argument
    if(present(i)) then
       i = hunt_list(x,xx,i)
       i1 = i
    else
       i1 = hunt_list(x,xx,0)
    end if

    if(i1 == 0) then
       i1 = 1
    else if(i1 == n) then
       i1 = n-1
    end if


    ! compute some parameters
    i2   = i+1
    x1   = x(i1)
    x2   = x(i2)
    f1   = f(i1)
    f2   = f(i2)
    fpp1 = fpp(i1)
    fpp2 = fpp(i2)
    h    = x2-x1
    a    = (x2-xx)/h
    b    = (xx-x1)/h

    if(present(ff)) then
       ! compute the function
       ff = a*f1 + b*f2 + ((a*a*a-a)*fpp1 +  (b*b*b-b)*fpp2)*(h*h)/6.0_dp
    endif

    if(present(ffp)) then
       ! compute the first derivative
       ffp = (f2-f1)/h  - (3.0_dp*a*a-1)*h*fpp1/6.0_dp  &
                        + (3.0_dp*b*b-1)*h*fpp2/6.0_dp
    end if

    if(present(ffpp)) then
       ! compute the second derivative
       ffpp = a*fpp1 + b*fpp2
    end if
    
    return
  end subroutine evaluate_cubic_spline


  
  subroutine parameters_cubic_spline(self,fpl,fpr)
    class(interp_1D_cubic), intent(inout) :: self
    real(dp), intent(in), optional :: fpl,fpr
    self%lnat = .not.present(fpl)
    self%rnat = .not.present(fpr)
    if(.not.self%lnat) self%fpl = fpl
    if(.not.self%rnat) self%fpr = fpr
    return
  end subroutine parameters_cubic_spline
  

  subroutine set_interp_1D_cubic(self,xx,ff)

    ! sets up a cubic spline from the data using natural
    ! boundary conditions.
    
    use module_constants
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), dimension(:), intent(in) :: xx,ff

    integer(i4b) :: i,n,info
    real(dp), dimension(:), allocatable :: dl,d,du
    real(dp), dimension(:,:), allocatable :: b
     
    ! set up the base type
    call self%interp_1D%set(xx,ff)
    
    ! dimension of the arrays
    n = self%n
    allocate(self%ff2(n))

    if(self%lnat .and. self%rnat) then
       call build_cubic_spline(n,self%xx,self%ff,self%ff2)
    end if
    if(self%lnat .and. .not.self%rnat) then
       call build_cubic_spline(n,self%xx,self%ff,self%ff2,fpr = self%fpr)
    end if
    if(.not.self%lnat .and. self%rnat) then
       call build_cubic_spline(n,self%xx,self%ff,self%ff2,fpl = self%fpl)
    end if
    if(.not.self%lnat .and. .not.self%rnat) then
       call build_cubic_spline(n,self%xx,self%ff,self%ff2,fpl = self%fpl,fpr = self%fpr)
    end if
    
    return
  end subroutine set_interp_1D_cubic
    

    
  function value_interp_1D_cubic(self,x) result(f)

    ! evaluates the cubic spline at a given location.
    ! Also computes and stores first and second derivatives, 
    ! with the results stored as type data  
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), intent(in) :: x

    integer(i4b) :: i1,i2
    real(dp) :: a,b,h,x1,x2,f1,f2,f21,f22,f

    
    ! get the indices and points
    call self%find(x,i1)
    i2 = i1+1

    ! compute some parameters
    x1  = self%xx(i1)
    x2  = self%xx(i2)
    f1  = self%ff(i1)
    f2  = self%ff(i2)
    f21 = self%ff2(i1)
    f22 = self%ff2(i2)
    h   = x2-x1
    a   = (x2-x)/h
    b   = (x-x1)/h

    ! compute the function  
    f = a*f1 + b*f2 + ((a*a*a-a)*f21 +  (b*b*b-b)*f22)*(h*h)/6.0_dp

    ! compute the first and second derivatives and store values   
    self%fp_save  = (f2-f1)/h  - (3.0_dp*a*a-1)*h*f21/6.0_dp  &
                          + (3.0_dp*b*b-1)*h*f22/6.0_dp
    self%fpp_save = a*f21 + b*f22
     
    
    return
  endfunction value_interp_1D_cubic


  function derivative_interp_1D_cubic(self,x) result(fp)

    ! evaluates the cubic spline at a given location.
    ! Also computes and stores first and second derivatives, 
    ! with the results stored as type data  
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), intent(in) :: x
    real(dp) :: f,fp

    if(x /= self%x_save) f = self%f(x)
    fp = self%fp_save

    return
  end function derivative_interp_1D_cubic


  function derivative2_interp_1D_cubic(self,x) result(fpp)

    ! evaluates the cubic spline at a given location.
    ! Also computes and stores first and second derivatives, 
    ! with the results stored as type data  
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), intent(in) :: x
    real(dp) :: f,fpp

    if(x /= self%x_save) f = self%f(x)
    fpp = self%fpp_save

    return
  end function derivative2_interp_1D_cubic
  
  
  subroutine deallocate_interp_1D_cubic(self)
    implicit none
    class(interp_1D_cubic) :: self
    deallocate(self%xx,self%ff,self%ff2)
    self%allocated = .false.    
    return
  end subroutine deallocate_interp_1D_cubic



  !==================================================================!
  !               routines for 2D bilinear interpolation             !
  !==================================================================!

  subroutine set_interp_2D(self,xx,yy,ff)
    implicit none
    class(interp_2D) :: self
    real(dp), dimension(:), intent(in) :: xx,yy
    real(dp), dimension(:,:), intent(in) :: ff

    
    integer(i4b) :: nx,ny,i,j
    
    ! check the size of the arrays
    nx = size(xx)
    ny = size(yy)
    call check(nx == size(ff,1) .and. ny == size(ff,2), &
               'set_linear_interp_2D','dimensions of arrays do not match')

    
    ! check the arrays are consecutive
    j = 0
    do i = 1,nx-1
       if(xx(i+1) >  xx(i)) then
          j = j+1
       else if(xx(i+1) > xx(i)) then
          j = j-1
       end if
    end do
    call check(j == nx-1 .or. j == -nx+1,'set_linear_interp_2D','x-array not consecutive')


    ! check the arrays are consecutive
    j = 0
    do i = 1,ny-1
       if(yy(i+1) >  yy(i)) then
          j = j+1
       else if(yy(i+1) > yy(i)) then
          j = j-1
       end if
    end do
    call check(j == ny-1 .or. j == -ny+1,'set_linear_interp_2D','y-array not consecutive')

    
    ! store the dimensions
    self%nx = nx
    self%ny = ny

    ! store the data
    self%xx = xx
    self%yy = yy
    self%ff = ff

    ! set dx_save
    self%dx_save = abs(xx(nx)-xx(1))/5.0_dp
    self%dy_save = abs(yy(ny)-yy(1))/5.0_dp

    ! build done
    self%allocated = .true.
    
    return
  end subroutine set_interp_2D


 
    
  subroutine deallocate_interp_2D(self)
    implicit none
    class(interp_2D) :: self
    deallocate(self%xx,self%yy,self%ff)
    self%allocated = .false.    
    return
  end subroutine deallocate_interp_2D


  subroutine find_interp_2D(self,x,y,ix,iy)
    use module_util, only: bisect_list,hunt_list
    implicit none
    class(interp_2D) :: self
    real(dp), intent(in) :: x,y
    integer(i4b), intent(out) :: ix,iy

    integer(i4b) :: il,iu,im,n

    ! is it worth hunting? 
    if(self%ix_save /= 0 .and. abs(x-self%x_save) > self%dx_save) then
       self%ix_save = 0
    end if

    ! find the point in x
    ix = hunt_list(self%xx,x,self%ix_save)

    ! is it worth hunting? 
    if(self%iy_save /= 0 .and. abs(y-self%y_save) > self%dy_save) then
       self%iy_save = 0
    end if

    ! find the point in y
    iy = hunt_list(self%yy,y,self%iy_save)

    ! _save the previous values
    self%ix_save = ix
    self%x_save = x
    self%iy_save = iy
    self%y_save = y
    
    return
  end subroutine find_interp_2D
  
  real(dp) function value_interp_2D(self,x,y) result(f)
    implicit none
    class(interp_2D) :: self
    real(dp), intent(in) :: x,y
    
    integer(i4b) :: ix,iy
    real(dp) :: x1,x2,y1,y2,f1,f2,f3,f4,s,t
    
    call self%find(x,y,ix,iy)    
    f =  bilinear_interp_calculation(self%xx,self%yy,self%ff,ix,iy,x,y) 
    
    return
  end function value_interp_2D


  real(dp) function bilinear_interp(xx,yy,ff,x,y,ix,iy) result(f)
    use module_util
    implicit none
    real(dp), dimension(:), intent(in) :: xx,yy
    real(dp), dimension(:,:), intent(in) :: ff
    real(dp), intent(in) :: x,y
    integer(i4b), intent(inout), optional :: ix,iy

    integer(i4b) :: jx,jy
    real(dp) :: x1,x2,y1,y2,f1,f2,f3,f4,s,t
    
    ! locate the indices
    if(present(ix) .and. present(iy)) then
       jx = hunt_list(xx,x,ix)
       jy = hunt_list(yy,y,iy)
       ix = jx
       iy = jy
    else
       jx = hunt_list(xx,x,0)
       jy = hunt_list(yy,y,0)       
    end if
    
    f =  bilinear_interp_calculation(xx,yy,ff,ix,iy,x,y) 
    
    return
  end function bilinear_interp


  real(dp) function bilinear_interp_calculation(xx,yy,ff,ix,iy,x,y) result(f)
    real(dp), dimension(:), intent(in) :: xx,yy
    real(dp), dimension(:,:), intent(in) :: ff
    integer(i4b), intent(in) :: ix,iy
    real(dp), intent(in) :: x,y
    real(dp) :: s,t,x1,x2,y1,y2,f1,f2,f3,f4

    x1 = xx(ix)
    x2 = xx(ix+1)
    y1 = yy(iy)
    y2 = yy(iy+1)

    f1 = ff(ix,iy)
    f2 = ff(ix+1,iy)
    f3 = ff(ix+1,iy+1)
    f4 = ff(ix,iy+1)
    
    s = (x-x1)/(x2-x1)
    t = (y-y1)/(y2-y1)
    f = (1-s)*(1-t)*f1 + s*(1-t)*f2 + s*t*f3 + (1-s)*t*f4
    
    return
  end function bilinear_interp_calculation


  !==================================================================!
  !                   routines for 2D bicubic spline                 !
  !==================================================================!


  subroutine set_interp_2D_bicubic_spline(self,xx,yy,ff)
    implicit none
    class(interp_2D_bicubic_spline) :: self
    real(dp), dimension(:), intent(in) :: xx,yy
    real(dp), dimension(:,:), intent(in) :: ff
    
    integer(i4b) :: i

    ! set up the basic type information
    call self%interp_2D%set(xx,yy,ff)

    ! build cubic splines along each column
    
    allocate(self%cs(self%ny))
    do i = 1,self%ny
       call self%cs(i)%set(xx,ff(:,i))	
    end do
    
    return
  end subroutine set_interp_2D_bicubic_spline


  subroutine deallocate_interp_2D_bicubic_spline(self)
    implicit none
    integer(i4b) :: i,n
    class(interp_2D_bicubic_spline) :: self
    deallocate(self%xx,self%yy,self%ff,self%cs)    
    self%allocated = .false.    
    return
  end subroutine deallocate_interp_2D_bicubic_spline


  real(dp) function value_interp_2D_bicubic_spline(self,x,y) result(f)
    implicit none
    class(interp_2D_bicubic_spline) :: self
    real(dp), intent(in) :: x,y
    
    integer(i4b) :: i
    real(dp), dimension(self%ny) :: g
    type(interp_1D_cubic) :: gs

    ! build temporary array for the given x
    do i = 1,self%ny
       g(i) = self%cs(i)%f(x)
    end do

    ! form the temporary cubic spline
    call gs%set(self%yy,g)

    ! evaluate the spline to get the result
    f = gs%f(y)

    
    return
  end function value_interp_2D_bicubic_spline



  
  
end module module_interp



