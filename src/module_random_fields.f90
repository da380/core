module module_random_fields

  use module_constants
  use module_util
  use module_error
  !use module_LAPACK
  use module_quadrature
  use module_interp
  use module_fftw3
  use module_SEM_1D
  use module_special_functions
  use module_spherical_harmonics
  use, intrinsic :: iso_c_binding
  implicit none


  !===================================================!
  !      types for random fields on an interval       !
  !===================================================!
  
  type, abstract :: GRF_1D
     integer(i4b) :: ns = 0
     real(dp) :: x1,x2
   contains
     procedure(GRF_1D_deallocate),    deferred :: deallocate
     procedure(GRF_1D_sample),    deferred :: sample
     procedure(GRF_1D_clear),     deferred :: clear
     procedure(GRF_1D_eval),      deferred :: eval
     procedure(GRF_1D_corr_eval), deferred :: corr_eval
  end type GRF_1D

  abstract interface
     
     subroutine GRF_1D_deallocate(self)
       import :: GRF_1D
       class(GRF_1D), intent(inout) :: self
     end subroutine GRF_1D_deallocate
     
     subroutine GRF_1D_sample(self,ns)
       use module_constants
       import :: GRF_1D
       class(GRF_1D), intent(inout) :: self
       integer(i4b), intent(in) :: ns
     end subroutine GRF_1D_sample
     
     subroutine GRF_1D_clear(self)
       import :: GRF_1D
       class(GRF_1D), intent(inout) :: self
     end subroutine GRF_1D_clear
     
     real(dp) function GRF_1D_eval(self,is,x)
       use module_constants
       import :: GRF_1D       
       class(GRF_1D), intent(inout) :: self
       integer(i4b), intent(in) :: is
       real(dp), intent(in) :: x
     end function GRF_1D_eval
     
     real(dp) function GRF_1D_corr_eval(self,x,c)
       use module_constants
       import :: GRF_1D       
       class(GRF_1D), intent(inout) :: self
       real(dp), intent(in) :: x
       real(dp), intent(in) :: c
     end function GRF_1D_corr_eval
     
  end interface

  type, extends(GRF_1D) :: GRF_1D_SEM
     logical :: allocated = .false.
     logical :: corr_point_set = .false.
     integer(i4b) :: ndim
     integer(i4b) :: mdim
     real(dp) :: c
     real(dp), dimension(:), allocatable :: Q
     real(dp), dimension(:), allocatable :: x
     real(dp), dimension(:,:), allocatable :: evec
     type(interp_1D_cubic), dimension(:), allocatable :: fun
     type(interp_1D_cubic) :: cfun
   contains
     procedure :: deallocate => deallocate_GRF_1D_SEM
     procedure :: sample => sample_GRF_1D_SEM
     procedure :: clear => clear_GRF_1D_SEM
     procedure :: eval => eval_GRF_1D_SEM
     procedure :: corr_eval => corr_eval_GRF_1D_SEM
  end type GRF_1D_SEM

  interface GRF_1D_SEM
     procedure :: build_GRF_1D_SEM
     procedure :: build_GRF_1D_radial_SEM
  end interface GRF_1D_SEM

  type, extends(GRF_1D) :: GRF_1D_Fourier
     logical :: allocated = .false.
     logical :: corr_point_set = .false.
     integer(i4b) :: n,i1,i2
     real(dp) :: c
     real(dp), dimension(:), allocatable :: x
     real(dp), dimension(:), allocatable :: Q
     type(C_PTR) :: plan_c2r
     type(interp_1D_cubic), dimension(:), allocatable :: fun
     type(interp_1D_cubic) :: cfun
   contains
     procedure :: deallocate => deallocate_GRF_1D_Fourier
     procedure :: sample => sample_GRF_1D_Fourier
     procedure :: clear => clear_GRF_1D_Fourier
     procedure :: eval => eval_GRF_1D_Fourier
     procedure :: corr_eval => corr_eval_GRF_1D_Fourier
  end type GRF_1D_Fourier

  interface GRF_1D_Fourier
          procedure :: build_GRF_1D_Fourier
  end interface GRF_1D_Fourier


  !===================================================!
  !            types for random fields on S2          !
  !===================================================!
  
  type, abstract :: GRF_S2
     integer(i4b) :: ns = 0
   contains
     procedure(GRF_S2_deallocate),    deferred :: deallocate
     procedure(GRF_S2_degree),    deferred :: degree
     procedure(GRF_S2_sample),    deferred :: sample
     procedure(GRF_S2_clear),     deferred :: clear
     procedure(GRF_S2_eval),      deferred :: eval
     procedure(GRF_S2_corr_eval), deferred :: corr_eval
     procedure(GRF_S2_coef),      deferred :: coef
     procedure(GRF_S2_corr_coef), deferred :: corr_coef
  end type GRF_S2

  abstract interface

     subroutine GRF_S2_deallocate(self)
       import :: GRF_S2
       class(GRF_S2), intent(inout) :: self
     end subroutine GRF_S2_deallocate
     
     integer(i4b) function GRF_S2_degree(self)
       use module_constants
       import :: GRF_S2
       class(GRF_S2), intent(inout) :: self
     end function GRF_S2_degree
     
     subroutine GRF_S2_sample(self,ns)
       use module_constants
       import :: GRF_S2       
       class(GRF_S2), intent(inout) :: self
       integer(i4b), intent(in) :: ns
     end subroutine GRF_S2_sample
     
     subroutine GRF_S2_clear(self)
       import :: GRF_S2       
       class(GRF_S2), intent(inout) :: self
     end subroutine GRF_S2_clear
     
     real(dp) function GRF_S2_eval(self,is,th,ph)
       use module_constants
       import :: GRF_S2
       class(GRF_S2), intent(inout) :: self
       integer(i4b), intent(in) :: is
       real(dp), intent(in) :: th,ph
     end function GRF_S2_eval
     
     real(dp) function GRF_S2_corr_eval(self,th,ph,th0,ph0)
       use module_constants
       import :: GRF_S2
       class(GRF_S2), intent(inout) :: self
       real(dp), intent(in) :: th,ph,th0,ph0
     end function GRF_S2_corr_eval
     
     subroutine GRF_S2_coef(self,is,lmax,ulm)
       use module_constants
       import :: GRF_S2
       class(GRF_S2), intent(inout) :: self
       integer(i4b), intent(in) :: is,lmax
       complex(dpc), dimension((lmax+1)*(lmax+2)/2), intent(out) :: ulm       
     end subroutine GRF_S2_coef
     
     subroutine GRF_S2_corr_coef(self,th0,ph0,lmax,clm)
       use module_constants
       import :: GRF_S2
       class(GRF_S2), intent(inout) :: self
       real(dp), intent(in) :: th0,ph0
       integer(i4b), intent(in) :: lmax
       complex(dpc), dimension((lmax+1)*(lmax+2)/2), intent(out) :: clm
     end subroutine GRF_S2_corr_coef
     
  end interface

  type, extends(GRF_S2) :: GRF_S2_SH
     logical :: allocated = .false.
     integer(i4b) :: lmax
     integer(i4b) :: ncoef
     real(dp), dimension(:), allocatable :: Q
     complex(dpc), dimension(:,:), allocatable :: ulm
   contains
     procedure :: deallocate    => deallocate_GRF_S2_SH
     procedure :: degree    => degree_GRF_S2_SH
     procedure :: sample    => sample_GRF_S2_SH
     procedure :: clear     => clear_GRF_S2_SH
     procedure :: eval      => eval_GRF_S2_SH
     procedure :: corr_eval => corr_eval_GRF_S2_SH
     procedure :: coef      => coef_GRF_S2_SH
     procedure :: corr_coef => corr_coef_GRF_S2_SH
  end type GRF_S2_SH
  

  !===================================================!
  !  types for random fields on a spherical annulus   !
  !===================================================!

  type, abstract :: GRF_B3
     integer(i4b) :: ns
     real(dp) :: a,b
   contains
     procedure(GRF_B3_deallocate),    deferred :: deallocate
     procedure(GRF_B3_degree),    deferred :: degree
     procedure(GRF_B3_sample),    deferred :: sample
     procedure(GRF_B3_clear),     deferred :: clear
     procedure(GRF_B3_eval),      deferred :: eval
     procedure(GRF_B3_corr_eval), deferred :: corr_eval
     procedure(GRF_B3_coef),      deferred :: coef
     procedure(GRF_B3_corr_coef), deferred :: corr_coef
  end type GRF_B3
  
  abstract interface
     
     subroutine GRF_B3_deallocate(self)
       import :: GRF_B3
       class(GRF_B3), intent(inout) :: self
     end subroutine GRF_B3_deallocate
     
     integer(i4b) function GRF_B3_degree(self)
       use module_constants
       import :: GRF_B3
       class(GRF_B3), intent(inout) :: self
     end function GRF_B3_degree
     
     subroutine GRF_B3_sample(self,ns)
       use module_constants
       import :: GRF_B3       
       class(GRF_B3), intent(inout) :: self
       integer(i4b), intent(in) :: ns
     end subroutine GRF_B3_sample
     
     subroutine GRF_B3_clear(self)
       import :: GRF_B3       
       class(GRF_B3), intent(inout) :: self
     end subroutine GRF_B3_clear
     
     real(dp) function GRF_B3_eval(self,is,r,th,ph)
       use module_constants
       import :: GRF_B3
       class(GRF_B3), intent(inout) :: self
       integer(i4b), intent(in) :: is
       real(dp), intent(in) :: r,th,ph
     end function GRF_B3_eval
     
     real(dp) function GRF_B3_corr_eval(self,r,th,ph,r0,th0,ph0)
       use module_constants
       import :: GRF_B3
       class(GRF_B3), intent(inout) :: self
       real(dp), intent(in) :: r,th,ph,r0,th0,ph0
     end function GRF_B3_corr_eval
     
     subroutine GRF_B3_coef(self,is,lmax,r,ulm)
       use module_constants
       import :: GRF_B3
       class(GRF_B3), intent(inout) :: self
       integer(i4b), intent(in) :: is,lmax
       real(dp), intent(in) :: r
       complex(dpc), dimension((lmax+1)*(lmax+2)/2), intent(out) :: ulm
     end subroutine GRF_B3_coef
     
     subroutine GRF_B3_corr_coef(self,r,r0,th0,ph0,lmax,ulm)
       use module_constants
       import :: GRF_B3
       class(GRF_B3), intent(inout) :: self
       real(dp), intent(in) :: r,r0,th0,ph0
       integer(i4b), intent(in) :: lmax
       complex(dpc), dimension((lmax+1)*(lmax+2)/2), intent(out) :: ulm
     end subroutine GRF_B3_corr_coef
     
  end interface

  type, extends(GRF_B3) :: GRF_B3_SEM
     logical :: allocated = .false.
     logical :: corr_point_set = .false.
     integer(i4b) :: lmax
     integer(i4b) :: ncoef
     real(dp) :: r0,th0,ph0
     class(GRF_1D), dimension(:), allocatable :: fun_l
   contains
     procedure :: deallocate    => deallocate_GRF_B3_SEM
     procedure :: degree    => degree_GRF_B3_SEM
     procedure :: sample    => sample_GRF_B3_SEM
     procedure :: clear     => clear_GRF_B3_SEM
     procedure :: eval      => eval_GRF_B3_SEM
     procedure :: corr_eval => corr_eval_GRF_B3_SEM
     procedure :: coef      => coef_GRF_B3_SEM
     procedure :: corr_coef => corr_coef_GRF_B3_SEM
  end type GRF_B3_SEM

contains


  !====================================================!
  !                 1D fields using SEM                !
  !====================================================!
  
  subroutine deallocate_GRF_1D_SEM(self)
    class(GRF_1D_SEM), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%Q,   & 
               self%x,    &          
               self%evec)
    call self%clear()    
    call self%cfun%deallocate()
    self%corr_point_set = .false.
    self%allocated = .false.
    return
  end subroutine deallocate_GRF_1D_SEM


  subroutine sample_GRF_1D_SEM(self,ns)
    class(GRF_1D_SEM), intent(inout) :: self
    integer(i4b), intent(in) :: ns
    integer(i4b) :: i,is
    real(dp), dimension(self%ndim) :: un
    real(dp), dimension(self%mdim) :: u
    call check(ns > 0,'sample_GRF_1D_SEM','sample size must be positive')
    if(self%ns /= ns) then
       call self%clear()
       allocate(self%fun(ns))
       self%ns = ns
    end if
    call random_seed()
    do is = 1,self%ns
       call normal_random_variable(un)
       un = un*sqrt(self%Q)
       u = 0.0_dp
       do i = 1,self%ndim
          u(:) = u(:) + un(i)*self%evec(:,i)
       end do
       call self%fun(is)%set(self%x,u)
    end do
    return
  end subroutine sample_GRF_1D_SEM


  subroutine clear_GRF_1D_SEM(self)
    class(GRF_1D_SEM), intent(inout) :: self
    integer(i4b) :: is 
    if(self%ns == 0) return
    do is = 1,self%ns
       call self%fun(is)%deallocate()
    end do
    deallocate(self%fun)
    self%ns = 0
    return
  end subroutine clear_GRF_1D_SEM
  
  real(dp) function eval_GRF_1D_SEM(self,is,x) result(f)
    class(GRF_1D_SEM), intent(inout) :: self
    integer(i4b), intent(in) :: is
    real(dp), intent(in) :: x
    f = self%fun(is)%f(x)
    return
  end function eval_GRF_1D_SEM

  real(dp) function corr_eval_GRF_1D_SEM(self,x,c) result(f)
    class(GRF_1D_SEM), intent(inout) :: self
    real(dp), intent(in) :: x,c
    integer(i4b) :: ic,i
    real(dp) :: c1,c2,v1,v2,v
    real(dp), dimension(self%mdim) :: u
    if(.not.self%corr_point_set .or. c /= self%c) then
       ic = bisect_list(self%x,c)
       c1 = self%x(ic)
       c2 = self%x(ic+1)
       u = 0.0_dp
       do i = 1,self%ndim
          v1 = self%evec(ic,i)
          v2 = self%evec(ic+1,i)
          v = v1 + (v2-v1)*(c-c1)/(c2-c1)
          u(:) = u(:) + self%Q(i)*v*self%evec(:,i)
       end do
       call self%cfun%set(self%x,u)       
       self%c = c
       self%corr_point_set = .true.
    end if
    f = self%cfun%f(x)
    return
  end function corr_eval_GRF_1D_SEM

  
  type(GRF_1D_SEM) function build_GRF_1D_SEM(x1,x2,lam,s,sig,eps) &
                            result(fun)
   
    real(dp), intent(in) :: x1,x2,lam,s,sig
    real(dp), intent(in), optional :: eps

    integer(i4b), parameter :: ngll = 5
    integer(i4b), parameter :: npad = 10
    real(dp), parameter :: eps_default = 1.e-5_dp

    integer(i4b) :: inode,ispec,count,kda,ldab,kdb,ldbb, &
                    ndim,i,j,k,jnode,knode,info,nmax,    &
                    nspec,ispec1,ispec2,i1,i2
    integer(i4b), dimension(:,:), allocatable :: ibool
    real(dp) :: x11,x22,xl,xr,dx,fac,sum,x,eps_loc,c
    real(dp), dimension(:), allocatable :: eval,work,jac
    real(dp), dimension(:,:), allocatable :: aa,bb,evec,hp,xx
    type(gauss_lobatto_quadrature) :: quad

        
    if(present(eps)) then
       eps_loc = eps
    else
       eps_loc = eps_default
    end if

    fun%x1 = x1
    fun%x2 = x2
    
    ! estimate the cutoff eigenvalue
    nmax = 0
    sum = 0.0_dp
    do       
       fac = 1.0_dp + (lam*pi*nmax/(x2-x1))**2
       fac = fac**(-s)
       sum = sum + fac
       if(fac/sum < eps_loc) exit
       nmax = nmax+1
    end do
    
    ! work out the mesh size
    dx = (x2-x1)/nmax
    nspec = (x2-x1)/dx
    dx = (x2-x1)/nspec
    nspec = nspec + 2*npad
    x11 = x1 - npad*dx
    x22 = x2 + npad*dx
    ispec1 = 1 + npad
    ispec2 = nspec - npad
   
   ! allocate the mesh arrays
    allocate(hp(ngll,ngll))
    allocate(jac(nspec))
    allocate(xx(ngll,nspec))
    allocate(ibool(ngll,nspec))

    
    ! get the gll points and weights  
    call quad%set(ngll)
    do inode = 1,ngll
       call lagrange_polynomial(quad%x(inode),ngll,quad%x,xx(:,1),hp(inode,:))
    end do
    
    ! build up the mesh
    xl = x11
    count = 0
    do ispec = 1,nspec
       xr = xl + dx
       jac(ispec) = 0.5_dp*(xr-xl)
       do inode = 1,ngll
          xx(inode,ispec) = xl + 0.5_dp*(xr-xl)*(quad%x(inode)+1.0_dp)
          count = count + 1
          ibool(inode,ispec) = count
       end do       
       xl = xr
       count = count-1
    end do
    
    ! allocate the matrices
    ndim = ibool(ngll,nspec)
    kda  = ngll-1
    ldab = kda+1
    kdb  = 0
    ldbb = kdb+1
    allocate(aa(ldab,ndim),bb(ldbb,ndim))
    aa = 0.0_dp
    bb = 0.0_dp
    
    ! build the matrices
    do ispec = 1,nspec
       do inode = 1,ngll
          i = ibool(inode,ispec)
          k = kdb+1
          bb(k,i) = bb(k,i) + quad%w(inode)*jac(ispec)
          do jnode = inode,ngll
             j = ibool(jnode,ispec)
             k = kda+1+i-j
             do knode = 1,ngll
                aa(k,j) = aa(k,j) + hp(knode,inode) &
                                  * hp(knode,jnode) &
                                  * quad%w(knode)   & 
                                  / jac(ispec)
             end do
          end do
       end do
    end do
    
    ! solve the eigenvalue problem
    allocate(eval(ndim))
    allocate(evec(ndim,ndim))
    allocate(work(3*ndim-2))
    call dsbgv('V','U',ndim,kda,kdb,aa,ldab,bb,ldbb,eval,evec,ndim,work,info)
    call check(info == 0,'set_gaussian_random_scalar_field_1D','problem with eigendecomposition')
    eval(1) = 0.0_dp

    
    ! work out eigenvalue cutoff
    fun%ndim = min(nmax,ndim)
    
    ! store the eigenvectors
    i1 = ibool(1,ispec1)
    i2 = ibool(ngll,ispec2)
    fun%mdim = i2-i1+1
   
    allocate(fun%evec(fun%mdim,fun%ndim),fun%x(fun%mdim))
    do i = 1,fun%ndim
       fun%evec(:,i) = evec(i1:i2,i)       
    end do
    do ispec = ispec1,ispec2
       do inode = 1,ngll
          i = ibool(inode,ispec)-i1+1
          fun%x(i) = xx(inode,ispec)
       end do
    end do
    
    ! set the covariance
    allocate(fun%Q(fun%ndim))
    fun%Q(:) = (1.0_dp + lam*lam*eval(1:fun%ndim))**(-s)
    j = fun%mdim/2
    sum = 0.0_dp
    do i = 1,fun%ndim
       sum = sum + fun%Q(i)*fun%evec(j,i)**2
    end do
    fun%Q = sig*sig*fun%Q/sum

    ! finish up
    fun%allocated = .true.
 
    return
  end function build_GRF_1D_SEM


  type(GRF_1D_SEM) function build_GRF_1D_radial_SEM(x1,x2,l,lam,s,sig,eps) &
                            result(fun)

    integer(i4b), intent(in) :: l
    real(dp), intent(in) :: x1,x2,lam,s,sig
    real(dp), intent(in), optional :: eps
    integer(i4b), parameter :: ngll = 5
    integer(i4b), parameter :: npad = 0
    real(dp), parameter :: eps_default = 1.e-5_dp

    integer(i4b) :: inode,ispec,count,kda,ldab,kdb,ldbb, &
                    ndim,i,j,k,jnode,knode,info,nmax,    &
                    nspec,ispec1,ispec2,i1,i2
    integer(i4b), dimension(:,:), allocatable :: ibool
    real(dp) :: x11,x22,xl,xr,dx,fac,sum,x,eps_loc,c
    real(dp), dimension(:), allocatable :: eval,work,jac
    real(dp), dimension(:,:), allocatable :: aa,bb,evec,hp,xx
    type(gauss_lobatto_quadrature) :: quad
        
    if(present(eps)) then
       eps_loc = eps
    else
       eps_loc = eps_default
   end if
   fun%x1 = x1
   fun%x2 = x2
   
    ! estimate the cutoff eigenvalue
    nmax = 0
   sum = 0.0_dp
    do       
       fac = 1.0_dp + (lam*pi*nmax/(x2-x1))**2
       fac = fac**(-s)
       sum = sum + fac
       if(fac/sum < eps_loc) exit
       nmax = nmax+1
    end do
    
    ! work out the mesh size
    dx = (x2-x1)/nmax
    nspec = (x2-x1)/dx
    dx = (x2-x1)/nspec
    nspec = nspec + 2*npad
    x11 = x1-npad*dx
    x22 = x2 + npad*dx
    ispec1 = 1 + npad
    ispec2 = nspec - npad
   
    ! allocate the mesh arrays
    allocate(hp(ngll,ngll))
    allocate(jac(nspec))
    allocate(xx(ngll,nspec))
    allocate(ibool(ngll,nspec))
    
    ! get the gll points and weights  
    call quad%set(ngll)
    do inode = 1,ngll
       call lagrange_polynomial(quad%x(inode),ngll,quad%x,xx(:,1),hp(inode,:))
    end do
    
    ! build up the mesh
    xl = x11
    count = 0
    do ispec = 1,nspec
       xr = xl + dx
       jac(ispec) = 0.5_dp*(xr-xl)
       do inode = 1,ngll
          xx(inode,ispec) = xl + 0.5_dp*(xr-xl)*(quad%x(inode)+1.0_dp)
          count = count + 1
          ibool(inode,ispec) = count
       end do       
       xl = xr
       count = count-1
    end do
    
    ! allocate the matrices
    ndim = ibool(ngll,nspec)
    kda  = ngll-1
    ldab = kda+1
    kdb  = 0
    ldbb = kdb+1
    allocate(aa(ldab,ndim),bb(ldbb,ndim))
    aa = 0.0_dp
    bb = 0.0_dp
    
    ! build the matrices
    do ispec = 1,nspec
       do inode = 1,ngll
          i = ibool(inode,ispec)
          k = kdb+1
          x = xx(inode,ispec)
          bb(k,i) = bb(k,i) + x*x*quad%w(inode)*jac(ispec)
          aa(k,i) = aa(k,i) + l*(l+1)*quad%w(inode)*jac(ispec)
          do jnode = inode,ngll
             j = ibool(jnode,ispec)
             k = kda+1+i-j
             do knode = 1,ngll
                x = xx(knode,ispec)
                aa(k,j) = aa(k,j) + x*x*hp(knode,inode) &
                                  * hp(knode,jnode) &
                                  * quad%w(knode)   & 
                                  / jac(ispec)
             end do
          end do
       end do
    end do
    
    ! solve the eigenvalue problem
    allocate(eval(ndim))
    allocate(evec(ndim,ndim))
    allocate(work(3*ndim-2))
    call dsbgv('V','U',ndim,kda,kdb,aa,ldab,bb,ldbb,eval,evec,ndim,work,info)
    call check(info == 0,'set_gaussian_random_scalar_field_1D','problem with eigendecomposition')
    eval(1) = 0.0_dp

    
    ! work out eigenvalue cutoff
    fun%ndim = min(nmax,ndim)
    
    ! store the eigenvectors
    i1 = ibool(1,ispec1)
    i2 = ibool(ngll,ispec2)
    fun%mdim = i2-i1+1
    
    allocate(fun%evec(fun%mdim,fun%ndim),fun%x(fun%mdim))
    do i = 1,fun%ndim
       fun%evec(:,i) = evec(i1:i2,i)       
    end do
    do ispec = ispec1,ispec2
       do inode = 1,ngll
          i = ibool(inode,ispec)-i1+1
          fun%x(i) = xx(inode,ispec)
       end do
    end do
    
    ! set the covariance
    allocate(fun%Q(fun%ndim))
    fun%Q(:) = (1.0_dp + lam*lam*eval(1:fun%ndim))**(-s)
    j = fun%mdim/2
    sum = 0.0_dp
    do i = 1,fun%ndim
       sum = sum + fun%Q(i)*fun%evec(j,i)**2
    end do
    fun%Q = sig*sig*fun%Q/sum

    ! finish up
    fun%allocated = .true.
 
    return
  end function build_GRF_1D_radial_SEM
  

  !====================================================!
  !                 1D fields using FFT                !
  !====================================================!
  
  
  subroutine deallocate_GRF_1D_Fourier(self)
    class(GRF_1D_Fourier), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%Q)
    call fftw_destroy_plan(self%plan_c2r)
    call self%clear()
    call self%cfun%deallocate()
    self%corr_point_set = .false.
    self%allocated = .false.
    return
  end subroutine deallocate_GRF_1D_Fourier


  subroutine sample_GRF_1D_Fourier(self,ns)
    class(GRF_1D_Fourier), intent(inout) :: self
    integer(i4b), intent(in) :: ns
    integer(i4b) :: i,n,i1,i2,is
    real(dp) :: r1,r2
    real(C_DOUBLE), pointer :: out(:)
    complex(C_DOUBLE_COMPLEX), pointer :: in(:)
    type(C_PTR) :: pin,pout
    call check(ns > 0,'sample_GRF_1D_Fourier','sample size must be positive')
    if(self%ns /= ns) then
       call self%clear()
       allocate(self%fun(ns))
       self%ns = ns
    end if    
    n = self%n
    i1 = self%i1
    i2 = self%i2
    pin = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    pout  = fftw_alloc_real(int(n, C_SIZE_T))
    call c_f_pointer(pin, in, [n/2+1])
    call c_f_pointer(pout,   out, [n])    
    call random_seed()
    do is = 1,self%ns
       call normal_random_variable(r1)
       in(1) = r1*sqrt(self%Q(1))
       do i = 2,n/2+1
          call normal_random_variable(r1,r2)
          in(i) = (r1+ii*r2)*sqrt(self%Q(i))
       end do
       call fftw_execute_dft_c2r(self%plan_c2r,in,out)
       call self%fun(is)%set(self%x(i1:i2),out(i1:i2))
    end do
    return
  end subroutine sample_GRF_1D_Fourier

  real(dp) function eval_GRF_1D_Fourier(self,is,x) result(f)
    class(GRF_1D_Fourier), intent(inout) :: self
    integer(i4b), intent(in) :: is
    real(dp), intent(in) :: x
    f = self%fun(is)%f(x)
    return
  end function eval_GRF_1D_Fourier  


  subroutine clear_GRF_1D_Fourier(self)
    class(GRF_1D_Fourier), intent(inout) :: self
    integer(i4b) :: is 
    if(self%ns == 0) return
    do is = 1,self%ns
       call self%fun(is)%deallocate()
    end do
    deallocate(self%fun)
    self%ns = 0
    return
  end subroutine clear_GRF_1D_Fourier

  
  real(dp) function corr_eval_GRF_1D_Fourier(self,x,c) result(f)
    class(GRF_1D_Fourier), intent(inout) :: self
    real(dp), intent(in) :: x,c

    integer(i4b) :: m,n,i1,i2
    real(dp) :: th,c1,c2
    real(C_DOUBLE), pointer :: out(:)
    complex(C_DOUBLE_COMPLEX), pointer :: in(:)
    type(C_PTR) :: pin,pout

    if(.not.self%corr_point_set .or. c /= self%c) then
       n = self%n
       i1 = self%i1
       i2 = self%i2
       c1 = self%x(1)
       c2 = self%x(n)
       th = twopi*(c-c1)/(c2-c1)
       pin = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
       pout  = fftw_alloc_real(int(n, C_SIZE_T))
       call c_f_pointer(pin, in, [n/2+1])
       call c_f_pointer(pout,   out, [n])
       in = 0.0_dp
       do m = 1,n/2+1
          in(m) = self%Q(m)*exp(-ii*(m-1)*th)
       end do
       call fftw_execute_dft_c2r(self%plan_c2r,in,out)
       call self%cfun%set(self%x(i1:i2),out(i1:i2))
       self%c = c
       self%corr_point_set = .true.
    end if
    f = self%cfun%f(x)
    return
  end function  corr_eval_GRF_1D_Fourier
  
  
  type(GRF_1D_Fourier) function build_GRF_1D_Fourier(x1,x2,lam,s, &
                                                    sig,eps,fftw_flag) result(fun)
    real(dp), intent(in) :: x1,x2,lam,s,sig
    real(dp), intent(in), optional :: eps
    integer(C_INT), intent(in), optional :: fftw_flag
    
    integer(i4b), parameter :: npad = 10
    real(dp), parameter :: eps_default = 1.0e-5_dp

    integer(i4b) :: i,n,i1,i2,ne
    real(dp) :: dx,x11,x22,k,eps_loc,fac,sum,dk

    real(C_DOUBLE), pointer :: out(:)
    complex(C_DOUBLE_COMPLEX), pointer :: in(:)
    type(C_PTR) :: pin,pout,plan
    integer(C_INT) :: plan_flag
    if(present(eps)) then
       eps_loc = eps
    else
       eps_loc = eps_default
    end if

    if(present(fftw_flag)) then
       plan_flag = fftw_flag
    else
       plan_flag = FFTW_MEASURE
    end if

    fun%x1 = x1
    fun%x2 = x2
    
    ! find maximum wavenumber
    i = 0
    dk = 1.0_dp/(10*lam)
    k = 0.0_dp
    sum = 0.0_dp    
    do
       fac = (1.0_dp+(lam*k)**2)**(-s)
       sum  = sum + fac
       if(fac/sum < eps_loc) exit
       k = k + dk
    end do
    dx = 0.5_dp/k

    ! number of points and range
    n = (x2-x1)/dx + 1
    dx = (x2-x1)/(n-1)
    n = n + 2*npad
    x11 = x1-npad*dx
    x22 = x2+npad*dx
    fun%i1 = 1+npad
    fun%i2 = n-npad

    ! extend to next power of 2
    ne = log(1.0_dp*n)/log(2.0_dp) + 1
    n = 2**ne
    fun%n = n

    ! set the points
    allocate(fun%x(n))
    do i = 1,n
       fun%x(i) = x11 + (i-1)*dx
    end do
    
    ! set the covariance
    allocate(fun%Q(n/2+1))
    sum = 0.0_dp
    do i = 1,n/2+1
       k = twopi*(i-1)/(dx*n)
       fac = (1.0_dp+(lam*k)**2)**(-s)
       if(i == 1) then
          sum = sum + fac
       else
          sum = sum + 2.0_dp*fac
       end if
       fun%Q(i) = fac
    end do
    fun%Q = sig*sig*fun%Q/sum

    ! set up fft plan
   pin = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    pout  = fftw_alloc_real(int(n, C_SIZE_T))
    call c_f_pointer(pin,   in, [n/2+1])
    call c_f_pointer(pout, out, [n])
    plan = fftw_plan_dft_c2r_1d(n,in,out,plan_flag)
    call fftw_free(pin)
    call fftw_free(pout)
    fun%plan_c2r = plan
    
    return
  end function build_GRF_1D_Fourier
  
  !====================================================!
  !              random fields on a sphere             !
  !====================================================!


  subroutine deallocate_GRF_S2_SH(self)
    class(GRF_S2_SH), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%Q)
    call self%clear()
    self%allocated = .false.
    return
  end subroutine deallocate_GRF_S2_SH


  integer(i4b) function degree_GRF_S2_SH(self) result(lmax)
    class(GRF_S2_SH), intent(inout) :: self
    lmax = self%lmax
    return
  end function degree_GRF_S2_SH
  
  subroutine sample_GRF_S2_SH(self,ns)
    class(GRF_S2_SH), intent(inout) :: self
    integer(i4b), intent(in) :: ns
    integer(i4b) :: l,m,ilm,is
    real(dp) :: r1,r2,std
    if(self%ns /= ns) then
       call self%clear()
       allocate(self%ulm(self%ncoef,ns))
       self%ns = ns
    end if
    call random_seed()
    do is = 1,self%ns
       self%ulm(:,is) = 0.0_dp
       ilm = 0
       do l = 0,self%lmax
          std = sqrt(self%Q(l))
          call normal_random_variable(r1)
          ilm = ilm + 1
          self%ulm(ilm,is) =  std*r1       
          do m = 1,l
             call normal_random_variable(r1,r2)
             ilm = ilm+1
             self%ulm(ilm,is) =  std*(r1+ii*r2)
          end do
       end do
       self%ulm(ilm,is) = 0.0_dp
    end do
    return
  end subroutine sample_GRF_S2_SH


  subroutine clear_GRF_S2_SH(self)
    class(GRF_S2_SH), intent(inout) :: self
    if(self%ns == 0) return
    deallocate(self%ulm)
    self%ns = 0
    return
  end subroutine clear_GRF_S2_SH
  
  real(dp) function eval_GRF_S2_SH(self,is,th,ph) result(u)
    class(GRF_S2_SH), intent(inout) :: self
    integer(i4b), intent(in) :: is
    real(dp), intent(in) :: th,ph
    integer(i4b) :: l,m,ilm
    real(dp) :: fac
    complex(dpc) :: ep,fp
    type(wigner_value) :: d
    u = 0.0_dp
    ilm = 0
    ep = exp(ii*ph)
    call d%init(th,0,self%lmax)
    do l = 0,self%lmax
       fac = sqrt((2*l+1)/fourpi)
       call d%next()
       ilm = ilm+1
       u  = u +  real(self%ulm(ilm,is)*fac*d%get(0,0),kind=dp)
       fp = 1.0_dp
       do m = 1,l
          fp = fp*ep
          ilm = ilm+1
          u = u + 2.0_dp*real(self%ulm(ilm,is)*fac*d%get(0,m)*fp,kind=dp)          
       end do
    end do        
    return
  end function eval_GRF_S2_SH

  
  real(dp) function corr_eval_GRF_S2_SH(self,th,ph,th0,ph0) result(u)
    class(GRF_S2_SH), intent(inout) :: self
    real(dp), intent(in) :: th,ph,th0,ph0
    integer(i4b) :: l
    real(dp) :: fac,del
    type(wigner_value) :: d
    del =  sin(th)*sin(th0)*cos(ph-ph0) + cos(th)*cos(th0) 
    del = acos(del)
    u = 0.0_dp
    call d%init(del,0,0)
    do l = 0,self%lmax
       fac = (2*l+1)/fourpi
       call d%next()
       u = u + fac*self%Q(l)*d%get(0,0)
    end do
    return
  end function corr_eval_GRF_S2_SH

  
  subroutine coef_GRF_S2_SH(self,is,lmax,ulm)
    class(GRF_S2_SH), intent(inout) :: self
    integer(i4b), intent(in) :: is,lmax
    complex(dpc), dimension((lmax+1)*(lmax+2)/2), intent(out) :: ulm
    integer(i4b) :: l,m,ilm
    ulm = 0.0_dp
    do l = 0,min(lmax,self%lmax)
       ilm = ilm + 1
       ulm(ilm) = self%ulm(ilm,is)
       do m = 1,l
          ilm = ilm+1
          ulm(ilm) = self%ulm(ilm,is)
       end do
    end do
    return
  end subroutine coef_GRF_S2_SH


  subroutine corr_coef_GRF_S2_SH(self,th0,ph0,lmax,clm)
    class(GRF_S2_SH), intent(inout) :: self
    real(dp), intent(in) :: th0,ph0
    integer(i4b), intent(in) :: lmax
    complex(dpc), dimension((lmax+1)*(lmax+2)/2), intent(out) :: clm
    integer(i4b) :: l,m,ilm,ll
    real(dp) :: fac
    complex(dpc) :: ep,fp
    type(wigner_value) :: d
    ll = min(lmax,self%lmax)
    clm = 0.0_dp
    ep = exp(-ii*ph0)
    call d%init(th0,0,ll)
    ilm = 0
    do l = 0,ll
       fac = sqrt((2*l+1)/fourpi)
       call d%next()
       ilm = ilm+1
       clm(ilm) = self%Q(l)*fac*d%get(0,0)
       fp = 1.0_dp
       do m = 1,l
          fp = fp*ep
          ilm = ilm+1
          clm(ilm) = self%Q(l)*fac*d%get(0,m)*fp
       end do       
    end do
    return
  end subroutine corr_coef_GRF_S2_SH

  
  
  type(GRF_S2_SH) function build_GRF_S2_SH(lam,s,sig,eps,asph) result(fun)
    real(dp), intent(in) :: lam
    real(dp), intent(in) :: s
    real(dp), intent(in) :: sig
    real(dp), intent(in), optional  :: eps
    logical, intent(in), optional :: asph 

    integer(i4b) :: l,l1,l2,ncoef
    real(dp), parameter :: eps_default = 1.e-5_dp
    real(dp) :: eps_loc,sum,fac
    
    
    if(present(eps)) then
       eps_loc = eps
    else
       eps_loc = eps_default
    end if

    ! work out lmax
    l1 = 0
    l2 = 2
    sum = 0.0_dp
    do 
       do l = l1,l2
          fac = (1.0_dp + lam*lam*l*(l+1))**(-s)          
          sum = sum + fac          
       end do
       if(fac/sum < eps_loc) exit
       l1 = l2+1
       l2 = 2*l2
    end do
    fun%lmax = l2
    fun%ncoef = (l2+1)*(l2+2)/2

    
    ! set the covariance
    allocate(fun%Q(0:fun%lmax))
    sum = 0.0_dp
    do l = 0,fun%lmax
       fun%Q(l) = (1.0_dp + lam*lam*l*(l+1))**(-s)
       sum = sum + (2*l+1)*fun%Q(l)/fourpi
    end do
    fun%Q = sig*sig*fun%Q/sum

    if(present(asph)) then
       if(asph) fun%Q(0) = 0.0_dp
    end if
   
    return
  end function build_GRF_S2_SH

  !====================================================!
  !        random fields in a spherical annulus        !
  !====================================================!


  subroutine deallocate_GRF_B3_SEM(self)
    class(GRF_B3_SEM), intent(inout) :: self
    integer(i4b) :: l
    call self%clear()
    self%corr_point_set = .false.
    self%allocated = .false.
    return
  end subroutine deallocate_GRF_B3_SEM

  integer(i4b) function degree_GRF_B3_SEM(self) result(lmax)
    class(GRF_B3_SEM), intent(inout) :: self
    lmax = self%lmax
    return
  end function degree_GRF_B3_SEM


  subroutine sample_GRF_B3_SEM(self,ns)
    class(GRF_B3_SEM), intent(inout) :: self
    integer(i4b), intent(in) :: ns
    integer(i4b) :: l
    do l = 0,self%lmax
       call self%fun_l(l)%sample(ns*(2*l+1))
    end do
    self%ns = ns
    return
  end subroutine sample_GRF_B3_SEM


  subroutine clear_GRF_B3_SEM(self)
    class(GRF_B3_SEM), intent(inout) :: self
    integer(i4b) :: l
    if(self%ns == 0) return
    do l = 0,self%lmax
       call self%fun_l(l)%clear()
    end do
    self%ns = 0
    return
  end subroutine clear_GRF_B3_SEM


  real(dp) function eval_GRF_B3_SEM(self,is,r,th,ph) result(u)
    class(GRF_B3_SEM), intent(inout) :: self
    real(dp), intent(in) :: r,th,ph
    integer(i4b), intent(in) :: is
    integer(i4b) :: l,m,ilm,jlm
    real(dp) :: fac
    complex(dpc) :: ep,fp,ulm
    type(wigner_value) :: d
    u = 0.0_dp
    ilm = 0
    ep = exp(ii*ph)
    call d%init(th,0,self%lmax)
    do l = 0,self%lmax
       fac = sqrt((2*l+1)/fourpi)
       call d%next()
       ilm = ilm+1
       jlm = (is-1)*(2*l+1)+1
       ulm = self%fun_l(l)%eval(jlm,r)
       u  = u +  real(ulm*fac*d%get(0,0),kind=dp)
       fp = 1.0_dp
       do m = 1,l
          fp = fp*ep
          ilm = ilm+1
          jlm = jlm+1
          ulm = self%fun_l(l)%eval(jlm,r)
          jlm = jlm+1
          ulm = ulm + ii*self%fun_l(l)%eval(jlm,r)
          u = u + 2.0_dp*real(ulm*fac*d%get(0,m)*fp,kind=dp)          
       end do
    end do            
    return
  end function eval_GRF_B3_SEM


  real(dp) function corr_eval_GRF_B3_SEM(self,r,th,ph,r0,th0,ph0) result(u)
    class(GRF_B3_SEM), intent(inout) :: self
    real(dp), intent(in) :: r,th,ph,r0,th0,ph0
    integer(i4b) :: l
    real(dp) :: del,fac,kl
    type(wigner_value) :: d
    del =  sin(th)*sin(th0)*cos(ph-ph0) + cos(th)*cos(th0) 
    del = acos(del)
    u = 0.0_dp
    call d%init(del,0,0)
    do l = 0,self%lmax
       fac = (2*l+1)/fourpi
       call d%next()
       kl = self%fun_l(l)%corr_eval(r,r0)
       u = u + fac*kl*d%get(0,0)
    end do    
    return
  end function corr_eval_GRF_B3_SEM

  subroutine coef_GRF_B3_SEM(self,is,lmax,r,ulm)
    class(GRF_B3_SEM), intent(inout) :: self
    integer(i4b), intent(in) :: is,lmax
    real(dp), intent(in) :: r
    complex(dpc), dimension((lmax+1)*(lmax+2)/2), intent(out) :: ulm
    integer(i4b) :: l,m,ilm,ll,jlm
    ll = min(lmax,self%lmax)
    ulm = 0.0_dp
    ilm = 0
    do l = 0,ll
       ilm = ilm+1
       jlm = (is-1)*(2*l+1) + 1
       ulm(ilm) = self%fun_l(l)%eval(jlm,r)
       do m = 1,l
          ilm = ilm + 1
          jlm = jlm+1
          ulm(ilm) = self%fun_l(l)%eval(jlm,r)
          jlm = jlm+1
          ulm(ilm) = ulm(ilm) + ii*self%fun_l(l)%eval(jlm,r)
       end do
    end do
    return
  end subroutine coef_GRF_B3_SEM


  subroutine corr_coef_GRF_B3_SEM(self,r,r0,th0,ph0,lmax,ulm)
    class(GRF_B3_SEM), intent(inout) :: self
    real(dp), intent(in) :: r,r0,th0,ph0
    integer(i4b), intent(in) :: lmax
    complex(dpc), dimension((lmax+1)*(lmax+2)/2), intent(out) :: ulm
    integer(i4b) :: l,m,ilm,ll
    real(dp) :: kl,fac
    complex(dpc) :: ep,fp
    type(wigner_value) :: d
    ll = min(lmax,self%lmax)
    ulm = 0.0_dp
    ilm = 0
    call d%init(th0,0,ll)
    ep = exp(-ii*ph0)
    do l = 0,ll
       ilm = ilm+1
       fac = sqrt((2*l+1)/fourpi)
       kl = self%fun_l(l)%corr_eval(r,r0)
       call d%next()
       ulm(ilm) = kl*fac*d%get(0,0)
       fp = 1.0_dp
       do m = 1,l
          ilm = ilm + 1
          fp = fp*ep
          ulm(ilm) = kl*fac*d%get(0,m)*fp
       end do
    end do
    return
  end subroutine corr_coef_GRF_B3_SEM
  
end module module_random_fields
