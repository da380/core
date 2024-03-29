module module_quadrature
  
  use module_constants
  implicit none
  
  type, abstract :: quadrature
     logical :: allocated = .false.
     integer(i4b) :: n
     real(dp) :: a,b
     real(dp), dimension(:), allocatable :: x
     real(dp), dimension(:), allocatable :: w
   contains
     procedure :: deallocate => deallocate_quadrature
     procedure :: allocate => allocate_quadrature   
     procedure :: order => order_quadrature     
     procedure :: points => points_quadrature
     procedure :: weights => weights_quadrature
     procedure :: trans => translate_quadrature
     procedure :: intr => int_quadrature_r
     procedure :: intc => int_quadrature_c
     generic :: int => intr,intc
  end type quadrature
  private :: int_quadrature_r,int_quadrature_c, &
             points_quadrature, weights_quadrature, &
             translate_quadrature,order_quadrature


  type, extends(quadrature) :: trapezoid_quadrature
   contains     
     procedure :: set => set_trapezoid_quadrature
  end type trapezoid_quadrature
  private :: set_trapezoid_quadrature
  
  type, extends(quadrature) :: gauss_quadrature
   contains     
     procedure :: set => set_gauss_quadrature
  end type gauss_quadrature
  private :: set_gauss_quadrature

  
  type, extends(quadrature) :: gauss_radau_quadrature
   contains     
     procedure :: set => set_gauss_radau_quadrature
  end type gauss_radau_quadrature
  private :: set_gauss_radau_quadrature


  type, extends(quadrature) :: gauss_lobatto_quadrature
   contains     
     procedure :: set => set_gauss_lobatto_quadrature
  end type gauss_lobatto_quadrature
  private :: set_gauss_lobatto_quadrature
  
  abstract interface
     function f_r_r(x) result(f)
       use module_constants
       implicit none
       real(dp), intent(in) :: x
       real(dp) :: f
     end function f_r_r     
     function f_r_c(x) result(f)
       use module_constants
       implicit none
       real(dp), intent(in) :: x
       complex(dpc) :: f
     end function f_r_c
  end interface
  
contains


  !=========================================================================!
  !                      general quadrature procedures                      !
  !=========================================================================!

  subroutine deallocate_quadrature(self)
    implicit none
    class(quadrature), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%x)
    deallocate(self%w)
    self%allocated = .false.
    return
  end subroutine deallocate_quadrature
  
  subroutine allocate_quadrature(self,n)
    use module_error
    implicit none
    class(quadrature), intent(inout) :: self
    integer(i4b), intent(in) :: n
    call check(n > 0,'allocate_quadrature','quadrature order less than 1')
    call self%deallocate()
    self%n = n
    allocate(self%x(n))
    allocate(self%w(n))
    self%allocated = .true.
    return
  end subroutine allocate_quadrature
  
  function order_quadrature(self) result(n)
    implicit none
    class(quadrature), intent(in) :: self
    integer(i4b) :: n
    n = self%n
    return
  end function order_quadrature

  function points_quadrature(self) result(x)
    implicit none
    class(quadrature), intent(in) :: self
    real(dp), dimension(:), allocatable :: x
    x = self%x
    return
  end function points_quadrature


  function weights_quadrature(self) result(w)
    implicit none
    class(quadrature), intent(in) :: self
    real(dp), dimension(:), allocatable :: w
    w = self%w
    return
  end function weights_quadrature
  
  function int_quadrature_r(self,f) result(int)
    implicit none
    class(quadrature), intent(inout) :: self
    procedure(f_r_r) :: f
    real(dp) :: int
    integer(i4b) :: i
    int = 0.0_dp
    do i = 1,self%n
       int = int + f(self%x(i))*self%w(i)
    end do    
    return
  end function int_quadrature_r

  function int_quadrature_c(self,f) result(int)
    implicit none
    class(quadrature), intent(inout) :: self
    procedure(f_r_c) :: f
    complex(dpc) :: int
    integer(i4b) :: i
    int = 0.0_dp
    do i = 1,self%n
       int = int + f(self%x(i))*self%w(i)
    end do    
    return
  end function int_quadrature_c
  
  subroutine translate_quadrature(self,c,d)
    implicit none
    class(quadrature), intent(inout) :: self
    real(dp), intent(in) :: c,d
    real(dp) :: jac,a,b
    a = self%a
    b = self%b
    jac = (d-c)/(b-a)
    self%x = (b-self%x)*c + (self%x-a)*d
    self%x = self%x/(b-a)
    self%w = self%w*jac
    self%a = c
    self%b = d        
    return
  end subroutine translate_quadrature


  !=========================================================================!
  !                      Specific quadrature procedures                     !
  !=========================================================================!


  subroutine set_trapezoid_quadrature(self,n,a,b)
    use module_error
    use module_special_functions
    implicit none    
    class(trapezoid_quadrature), intent(inout) :: self
    integer(i4b), intent(in) :: n
    real(dp), intent(in) :: a,b
    integer(i4b) :: i
    real(dp) :: dx
    call check(n > 1,'set_trapezoid_quadrature','order must be at least 2')
    call self%allocate(n)
    self%a = a
    self%b = b
    dx = (b-a)/(n-1)
    do i = 1,n
       self%x(i) = a+(i-1)*dx
    end do
    self%w = 1.0_dp
    self%w(1) = 0.5_dp
    self%w(n) = 0.5_dp
    return
  end subroutine set_trapezoid_quadrature
  
  
  subroutine set_gauss_quadrature(self,n)
    use module_special_functions
    implicit none    
    class(gauss_quadrature), intent(inout) :: self
    integer(i4b), intent(in) :: n
    call self%allocate(n)    
    self%a = -1.0_dp
    self%b =  1.0_dp
    call build_gauss_quadrature(n,self%x,self%w)
    return
  end subroutine set_gauss_quadrature

  subroutine build_gauss_quadrature(n,x,w)
    use module_special_functions
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), dimension(n), intent(out) :: x,w

    integer(i4b) :: il,iu,abstol,m,ldz,lwork,liwork,info,i
    integer(i4b), dimension(10*n) :: iwork
    integer(i4b), dimension(2*n) :: isuppz
    real(dp) :: vl,vu,ix
    real(dp), dimension(n) :: alpha,beta
    real(dp), dimension(18*n) :: work
    real(dp), dimension(n,n) :: z
    
    ! set up the matrix
    do i = 1,n
       ix = i
       alpha(i) = 0.0_dp
       beta(i)  = ix*ix/((2*ix-1)*(2*ix+1))
       beta(i)  = sqrt(beta(i))
    end do
    
    ! solve the eigenvalue problem
    lwork  = 18*n
    liwork = 10*n
    call dstegr('V','A',n,alpha,beta,vl,vu,il,iu,abstol, &
                 m,x,z,n,isuppz,work,lwork,iwork,liwork,info)		

    ! get the weights
    w = 2.0_dp*z(1,:)**2
    
    return
  end subroutine build_gauss_quadrature


  subroutine set_gauss_radau_quadrature(self,n,right)
    use module_error
    use module_special_functions
    implicit none    
    class(gauss_radau_quadrature), intent(inout) :: self
    integer(i4b), intent(in) :: n
    logical, intent(in), optional :: right
    call check(n > 1,'set_gauss_radau_quadrature','order must be at least 2')
    call self%allocate(n)
    self%a = -1.0_dp
    self%b =  1.0_dp    
    call build_gauss_radau_quadrature(n-1,self%x,self%w,right)
    return
  end subroutine set_gauss_radau_quadrature
  
  subroutine build_gauss_radau_quadrature(n,x,w,right)
    use module_special_functions
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), dimension(n+1), intent(out) :: x,w
    logical, intent(in), optional :: right

    logical :: do_right
    integer(i4b) :: il,iu,abstol,m,ldz,lwork,liwork,info,i
    integer(i4b), dimension(10*(n+1)) :: iwork
    integer(i4b), dimension(2*(n+1)) :: isuppz
    real(dp) :: vl,vu,ix,x1,x2
    real(dp), dimension(n) :: alpha,beta
    real(dp), dimension(n+1) :: d,e
    real(dp), dimension(n,1) :: b
    real(dp), dimension(18*(n+1)) :: work
    real(dp), dimension(n+1,n+1) :: z

    x1 = -1.0_dp
    x2 =  1.0_dp
    
    ! set up the matrix
    do i = 1,n
       ix = i
       alpha(i) = 0.0_dp
       beta(i)  = ix*ix/((2*ix-1)*(2*ix+1))
       beta(i)  = sqrt(beta(i))
    end do

    
    if(present(right)) then
       do_right = right
    else
       do_right = .false.
    end if

    ! set up the linear system to be solved
    if(do_right) then       
       d(1:n)   = x2-alpha(1:n)
       e(1:n-1) = -beta(1:n-1)
       b = 0.0_dp
       b(n,1) = -beta(n)**2
    else
       d(1:n)   = alpha(1:n)-x1
       e(1:n-1) = beta(1:n-1)
       b = 0.0_dp
       b(n,1) = beta(n)**2
    end if
    
    ! solve the linear system
    call  dptsv(n,1,d,e,b,n,info)

    ! set up the matrix for the eigenvalue problem
    d(1:n) = alpha
    if(do_right) then
       d(n+1) = x2 + b(n,1)
    else
       d(n+1) = x1 + b(n,1)
    end if
    e(1:n) = beta
    e(n+1) = 0.0_dp
    
    ! solve the eigenvalue problem
    lwork  = 18*(n+1)
    liwork = 10*(n+1)
    call dstegr('V','A',n+1,d,e,vl,vu,il,iu,abstol, &
                 m,x,z,n+1,isuppz,work,lwork,iwork,liwork,info)		

    ! get the weights
    w = 2.0_dp*z(1,:)**2

    if(do_right) then
       x(n+1) = 1.0_dp
    else
       x(1)  = -1.0_dp
    end if
    
    return
  end subroutine build_gauss_radau_quadrature


  subroutine set_gauss_lobatto_quadrature(self,n)
    use module_error
    use module_special_functions
    implicit none    
    class(gauss_lobatto_quadrature), intent(inout) :: self
    integer(i4b), intent(in) :: n
    call check(n > 1,'set_gauss_lobatto_quadrature','order must be at least 2')
    call self%allocate(n)
    self%a = -1.0_dp
    self%b =  1.0_dp  
    call build_gauss_lobatto_quadrature(n-1,self%x,self%w)
    return
  end subroutine set_gauss_lobatto_quadrature
  


  subroutine build_gauss_lobatto_quadrature(n,x,w)
    use module_special_functions
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), dimension(n+1), intent(out) :: x,w

    logical :: do_right
    integer(i4b) :: il,iu,abstol,m,ldz,lwork,liwork,info,i
    integer(i4b), dimension(10*(n+1)) :: iwork
    integer(i4b), dimension(2*(n+1)) :: isuppz
    real(dp) :: vl,vu,gamma,mu,x1,x2,mu0,ix
    real(dp), dimension(n) :: alpha,beta
    real(dp), dimension(n+1) :: d,e
    real(dp), dimension(n,1) :: b
    real(dp), dimension(18*(n+1)) :: work
    real(dp), dimension(n+1,n+1) :: z

    ! get end points
    x1 = -1.0_dp
    x2 =  1.0_dp

    ! first moment
    mu0 = 2.0_dp
    
    ! set up the matrix
    do i = 1,n
       ix = i
       alpha(i) = 0.0_dp
       beta(i)  = ix*ix/((2*ix-1)*(2*ix+1))
       beta(i)  = sqrt(beta(i))
    end do

    ! set up and solve the first linear system
    d(1:n)   = alpha(1:n)-x1
    e(1:n-1) = beta(1:n-1)
    b        = 0.0_dp
    b(n,1)   = 1.0_dp
    call  dptsv(n,1,d,e,b,n,info)
    gamma = b(n,1)

    ! set up and solve the second linear system
    d(1:n)   = x2-alpha(1:n)
    e(1:n-1) = -beta(1:n-1)
    b = 0.0_dp
    b(n,1) = -1.0_dp
    call  dptsv(n,1,d,e,b,n,info)
    mu = b(n,1)

    ! set up and solve the eigenvalue problem
    d(1:n) = alpha
    e(1:n-1) = beta(1:n-1)
    e(n)   = sqrt((x2-x1)/(gamma-mu))
    d(n+1) = x1 + gamma*e(n)**2      
    lwork  = 18*(n+1)
    liwork = 10*(n+1)
    call dstegr('V','A',n+1,d,e,vl,vu,il,iu,abstol, &
                 m,x,z,n+1,isuppz,work,lwork,iwork,liwork,info)		

    ! get the weights
    w = mu0*z(1,:)**2

    ! force end points to be exact
    x(1) = -1.0_dp
    x(n+1) = 1.0_dp
    
    return
  end subroutine build_gauss_lobatto_quadrature

  
end module module_quadrature
