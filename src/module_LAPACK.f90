module module_LAPACK

  use module_constants
  use module_error
  implicit none

  type rmat
     ! basic data
     logical :: allocated = .false.
     logical :: square = .false.
     logical :: check = .false.
     integer(i4b) :: m = 0
     integer(i4b) :: n = 0
     real(dp), dimension(:,:), allocatable :: elem
     ! LU data
     logical :: LU_factorised = .false.
     integer(i4b), dimension(:), allocatable :: ipiv
   contains
     procedure :: deallocate => deallocate_rmat
     procedure :: inc        => increment_rmat
     procedure :: LU         => LU_rmat
     procedure :: LU_backsub => LU_backsub_rmat
  end type rmat

  interface rmat
     procedure :: allocate_rmat_from_dimensions
     procedure :: allocate_square_rmat_from_dimensions
     procedure :: allocate_rmat_from_rmat
     procedure :: allocate_rmat_from_array
  end interface rmat

  type, extends(rmat) :: brmat
     integer(i4b) :: kl
     integer(i4b) :: ku
     integer(i4b) :: ld
   contains
     procedure :: deallocate => deallocate_brmat
     procedure :: inc        => increment_brmat
     procedure :: LU         => LU_brmat
     procedure :: LU_backsub => LU_backsub_brmat
  end type brmat


  interface brmat
     procedure :: allocate_brmat_from_dimensions
     procedure :: allocate_square_brmat_from_dimensions
     procedure :: allocate_brmat_from_brmat
  end interface brmat


  type, extends(rmat) :: sbrmat
     integer(i4b) :: kl
     integer(i4b) :: ku
     integer(i4b) :: ld
   contains
     procedure :: deallocate => deallocate_sbrmat
     procedure :: inc        => increment_sbrmat
     procedure :: LU         => LU_sbrmat
     procedure :: LU_backsub => LU_backsub_sbrmat
  end type sbrmat


  interface sbrmat
     procedure :: allocate_sbrmat_from_dimensions
     procedure :: allocate_square_sbrmat_from_dimensions
     procedure :: allocate_sbrmat_from_sbrmat
  end interface sbrmat



  type complex_mat
     ! basic data
     logical :: allocated = .false.
     logical :: square = .false.
     logical :: check = .false.
     integer(i4b) :: m = 0
     integer(i4b) :: n = 0
     complex(dpc), dimension(:,:), allocatable :: elem
     ! LU data
     logical :: LU_factorised = .false.
     integer(i4b), dimension(:), allocatable :: ipiv
   contains
     procedure :: deallocate => deallocate_complex_mat
     procedure :: increment_real_complex_mat
     procedure :: increment_complex_complex_mat
     generic   :: inc => increment_real_complex_mat, &
                         increment_complex_complex_mat
     procedure :: LU         => LU_complex_mat
     procedure :: LU_backsub => LU_backsub_complex_mat
  end type complex_mat

  interface complex_mat
     procedure :: allocate_complex_mat_from_dimensions
     procedure :: allocate_square_complex_mat_from_dimensions
     procedure :: allocate_complex_mat_from_complex_mat
     procedure :: allocate_complex_mat_from_rmat
     procedure :: allocate_complex_mat_from_array
  end interface complex_mat
  
contains


  !===============================================!
  !          procedures for real matrices         !
  !===============================================!

  !---------------------------------------!
  !       methods for real matrices       !
  !---------------------------------------!
    
  subroutine deallocate_rmat(self)
    class(rmat), intent(inout) :: self
    if(.not.self%allocated) return
    self%m = 0
    self%n = 0
    deallocate(self%elem)
    if(self%LU_factorised) then
       deallocate(self%ipiv)
       self%LU_factorised = .false.
    end if
    self%allocated = .false.    
    return
  end subroutine deallocate_rmat
  
  subroutine increment_rmat(self,i,j,val)
    class(rmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'increment_rmat','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'increment_rmat','row index out of range')
       call check(j >= 1 .and. j <= self%n,'increment_rmat','column index out of range')
    end if
    self%elem(i,j) = self%elem(i,j) + val
    return
  end subroutine increment_rmat
    
  subroutine LU_rmat(self)
    class(rmat), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'LU_rmat','matrix not allocated')
    if(self%LU_factorised) return
    allocate(self%ipiv(min(self%m,self%n)))
    call dgetrf(self%m,self%n,self%elem,self%m,self%ipiv,info)
    call check(info == 0,'LU_rmat','problem with factorisation')
    self%LU_factorised = .true.
    return
  end subroutine LU_rmat

  subroutine LU_backsub_rmat(self,b,trans)
    class(rmat), intent(in) :: self
    class(rmat), intent(inout) :: b
    logical, intent(in), optional :: trans
    character(len=1) :: trans_char
    integer(i4b) :: m,nrhs,info
    trans_char = 'N'    
    if(present(trans)) then
       if(trans) trans_char = 'T'
    end if
    m = self%m
    nrhs = b%n
    call check(self%square,'LU_backsub_rmat','matrix not square')
    call check(self%LU_factorised,'LU_backsub_rmat','matrix not factorised')
    call check(m == b%m,'LU_backsub_rmat','rhs has the wrong dimensions')
    call dgetrs	(trans_char,m,nrhs,self%elem,m,self%ipiv,b%elem,m,info)
    call check(info == 0,'LU_backsub_rmat','problem with back substitution')
    return
  end subroutine LU_backsub_rmat


  
  !---------------------------------------!
  !     constructors for real matrices    !
  !---------------------------------------!

  type(rmat) function allocate_rmat_from_dimensions(m,n) result(a)
    integer(i4b), intent(in) :: m,n
    a%m = m
    a%n = n
    allocate(a%elem(m,n))
    a%elem = 0.0_dp
    if(m == n) a%square = .true.
    a%allocated = .true.
    return
  end function allocate_rmat_from_dimensions


  type(rmat) function allocate_square_rmat_from_dimensions(m) result(a)
    integer(i4b), intent(in) :: m    
    a = rmat(m,m)
    return
  end function allocate_square_rmat_from_dimensions


  type(rmat) function allocate_rmat_from_rmat(b) result(a)
    type(rmat), intent(in) :: b
    a%m = b%m
    a%n = b%n
    allocate(a%elem,source = b%elem)
    a%square = b%square
    a%allocated = .true.
    return
  end function allocate_rmat_from_rmat


  type(rmat) function allocate_rmat_from_array(b) result(a)
    real(dp), dimension(:,:), intent(in) :: b
    a%m = size(b,1)
    a%n = size(b,2)
    if(a%m == a%n) a%square = .true.
    allocate(a%elem,source = b)
    a%allocated = .true.
    return
  end function allocate_rmat_from_array


  type(rmat) function unit_rmat(m) result(a)
    integer(i4b), intent(in) :: m
    integer(i4b) :: i
    a = rmat(m)
    do i = 1,m
       call a%inc(i,i,1.0_dp)
    end do
    return
  end function unit_rmat


  !===============================================!
  !      procedures for band real matrices      !
  !===============================================!

  !---------------------------------------!
  !    methods for band real matrices   !
  !---------------------------------------!
    
  subroutine deallocate_brmat(self)
    class(brmat), intent(inout) :: self
    if(.not.self%allocated) return
    self%m = 0
    self%n = 0
    self%kl = 0
    self%ku = 0
    self%ld = 0
    deallocate(self%elem)
    if(self%LU_factorised) then
       deallocate(self%ipiv)
       self%LU_factorised = .false.
    end if
    self%allocated = .false.    
    return
  end subroutine deallocate_brmat
  
  subroutine increment_brmat(self,i,j,val)
    class(brmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'increment_brmat','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'increment_brmat','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'increment_brmat','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    self%elem(k,j) = self%elem(k,j) + val
    return
  end subroutine increment_brmat
    
  subroutine LU_brmat(self)
    class(brmat), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'LU_brmat','matrix not allocated')
    if(self%LU_factorised) return
    allocate(self%ipiv(min(self%m,self%n)))
    call dgbtrf(self%m,self%n,self%kl,self%ku,self%elem,self%ld,self%ipiv,info)
    call check(info == 0,'LU_brmat','problem with factorisation')
    self%LU_factorised = .true.
    return
  end subroutine LU_brmat

  subroutine LU_backsub_brmat(self,b,trans)
    class(brmat), intent(in) :: self
    class(rmat), intent(inout) :: b
    logical, intent(in), optional :: trans
    character(len=1) :: trans_char
    integer(i4b) :: m,nrhs,info
    trans_char = 'N'    
    if(present(trans)) then
       if(trans) trans_char = 'T'
    end if
    m = self%m
    nrhs = b%n
    call check(self%square,'LU_backsub_brmat','matrix not square')
    call check(self%LU_factorised,'LU_backsub_brmat','matrix not factorised')
    call check(m == b%m,'LU_backsub_brmat','rhs has the wrong dimensions')
    call dgbtrs(trans_char,m,self%kl,self%ku,nrhs,self%elem,self%ld,self%ipiv,b%elem,m,info)	
    call check(info == 0,'LU_backsub_brmat','problem with back substitution')
    return
  end subroutine LU_backsub_brmat


  !----------------------------------------------!
  !     constructors for band real matrices    !
  !----------------------------------------------!

  type(brmat) function allocate_brmat_from_dimensions(m,n,kl,ku) result(a)
    integer(i4b), intent(in) :: m,n,kl,ku
    a%m  = m
    a%n  = n
    a%kl = kl
    a%ku = ku
    a%ld = 2*a%kl + a%ku + 1
    allocate(a%elem(a%ld,n))
    a%elem = 0.0_dp
    if(m == n) a%square = .true.
    a%allocated = .true.
    return
  end function allocate_brmat_from_dimensions


  type(brmat) function allocate_square_brmat_from_dimensions(m,kl,ku) result(a)
    integer(i4b), intent(in) :: m,kl,ku    
    a = brmat(m,m,kl,ku)
    return
  end function allocate_square_brmat_from_dimensions


  type(brmat) function allocate_brmat_from_brmat(b) result(a)
    type(brmat), intent(in) :: b
    a%m = b%m
    a%n = b%n
    a%kl = b%kl
    a%ku = b%ku
    a%ld = b%ld
    allocate(a%elem,source = b%elem)
    a%square = b%square
    a%allocated = .true.
    return
  end function allocate_brmat_from_brmat

  
  type(brmat) function unit_brmat(m,kl,ku) result(a)
    integer(i4b), intent(in) :: m,kl,ku
    integer(i4b) :: i
    a = brmat(m,kl,ku)
    do i = 1,m
       call a%inc(i,i,1.0_dp)
    end do
    return
  end function unit_brmat



  !===============================================!
  !      procedures for symmetric_band real matrices      !
  !===============================================!

  !---------------------------------------!
  !    methods for symmetric_band real matrices   !
  !---------------------------------------!
    
  subroutine deallocate_sbrmat(self)
    class(sbrmat), intent(inout) :: self
    if(.not.self%allocated) return
    self%m = 0
    self%n = 0
    self%kl = 0
    self%ku = 0
    self%ld = 0
    deallocate(self%elem)
    if(self%LU_factorised) then
       deallocate(self%ipiv)
       self%LU_factorised = .false.
    end if
    self%allocated = .false.    
    return
  end subroutine deallocate_sbrmat
  
  subroutine increment_sbrmat(self,i,j,val)
    class(sbrmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'increment_sbrmat','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'increment_sbrmat','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'increment_sbrmat','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    self%elem(k,j) = self%elem(k,j) + val
    return
  end subroutine increment_sbrmat
    
  subroutine LU_sbrmat(self)
    class(sbrmat), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'LU_sbrmat','matrix not allocated')
    if(self%LU_factorised) return
    allocate(self%ipiv(min(self%m,self%n)))
    call dgbtrf(self%m,self%n,self%kl,self%ku,self%elem,self%ld,self%ipiv,info)
    call check(info == 0,'LU_sbrmat','problem with factorisation')
    self%LU_factorised = .true.
    return
  end subroutine LU_sbrmat

  subroutine LU_backsub_sbrmat(self,b,trans)
    class(sbrmat), intent(in) :: self
    class(rmat), intent(inout) :: b
    logical, intent(in), optional :: trans
    character(len=1) :: trans_char
    integer(i4b) :: m,nrhs,info
    trans_char = 'N'    
    if(present(trans)) then
       if(trans) trans_char = 'T'
    end if
    m = self%m
    nrhs = b%n
    call check(self%square,'LU_backsub_sbrmat','matrix not square')
    call check(self%LU_factorised,'LU_backsub_sbrmat','matrix not factorised')
    call check(m == b%m,'LU_backsub_sbrmat','rhs has the wrong dimensions')
    call dgbtrs(trans_char,m,self%kl,self%ku,nrhs,self%elem,self%ld,self%ipiv,b%elem,m,info)	
    call check(info == 0,'LU_backsub_sbrmat','problem with back substitution')
    return
  end subroutine LU_backsub_sbrmat


  !----------------------------------------------!
  !     constructors for symmetric_band real matrices    !
  !----------------------------------------------!

  type(sbrmat) function allocate_sbrmat_from_dimensions(m,n,kl,ku) result(a)
    integer(i4b), intent(in) :: m,n,kl,ku
    a%m  = m
    a%n  = n
    a%kl = kl
    a%ku = ku
    a%ld = 2*a%kl + a%ku + 1
    allocate(a%elem(a%ld,n))
    a%elem = 0.0_dp
    if(m == n) a%square = .true.
    a%allocated = .true.
    return
  end function allocate_sbrmat_from_dimensions


  type(sbrmat) function allocate_square_sbrmat_from_dimensions(m,kl,ku) result(a)
    integer(i4b), intent(in) :: m,kl,ku    
    a = sbrmat(m,m,kl,ku)
    return
  end function allocate_square_sbrmat_from_dimensions


  type(sbrmat) function allocate_sbrmat_from_sbrmat(b) result(a)
    type(sbrmat), intent(in) :: b
    a%m = b%m
    a%n = b%n
    a%kl = b%kl
    a%ku = b%ku
    a%ld = b%ld
    allocate(a%elem,source = b%elem)
    a%square = b%square
    a%allocated = .true.
    return
  end function allocate_sbrmat_from_sbrmat

  
  type(sbrmat) function unit_sbrmat(m,kl,ku) result(a)
    integer(i4b), intent(in) :: m,kl,ku
    integer(i4b) :: i
    a = sbrmat(m,kl,ku)
    do i = 1,m
       call a%inc(i,i,1.0_dp)
    end do
    return
  end function unit_sbrmat



  


  !===============================================!
  !        procedures for complex matrices        !
  !===============================================!

  !---------------------------------------!
  !      methods for complex matrices     !
  !---------------------------------------!
    
  subroutine deallocate_complex_mat(self)
    class(complex_mat), intent(inout) :: self
    if(.not.self%allocated) return
    self%m = 0
    self%n = 0
    deallocate(self%elem)
    if(self%LU_factorised) then
       deallocate(self%ipiv)
       self%LU_factorised = .false.
    end if
    self%allocated = .false.    
    return
  end subroutine deallocate_complex_mat
  
  subroutine increment_complex_complex_mat(self,i,j,val)
    class(complex_mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'increment_complex_mat','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'increment_complex_mat','row index out of range')
       call check(j >= 1 .and. j <= self%n,'increment_complex_mat','column index out of range')
    end if
    self%elem(i,j) = self%elem(i,j) + val
    return
  end subroutine increment_complex_complex_mat

  subroutine increment_real_complex_mat(self,i,j,val)
    class(complex_mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'increment_complex_mat','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'increment_complex_mat','row index out of range')
       call check(j >= 1 .and. j <= self%n,'increment_complex_mat','column index out of range')
    end if
    self%elem(i,j) = self%elem(i,j) + val
    return
  end subroutine increment_real_complex_mat
    
  subroutine LU_complex_mat(self)
    class(complex_mat), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'LU_complex_mat','matrix not allocated')
    if(self%LU_factorised) return
    allocate(self%ipiv(min(self%m,self%n)))
    call zgetrf(self%m,self%n,self%elem,self%m,self%ipiv,info)
    call check(info == 0,'LU_complex_mat','problem with factorisation')
    self%LU_factorised = .true.
    return
  end subroutine LU_complex_mat

  subroutine LU_backsub_complex_mat(self,b,trans)
    class(complex_mat), intent(in) :: self
    class(complex_mat), intent(inout) :: b
    logical, intent(in), optional :: trans
    character(len=1) :: trans_char
    integer(i4b) :: m,nrhs,info
    trans_char = 'N'    
    if(present(trans)) then
       if(trans) trans_char = 'T'
    end if
    m = self%m
    nrhs = b%n
    call check(self%square,'LU_backsub_complex_mat','matrix not square')
    call check(self%LU_factorised,'LU_backsub_complex_mat','matrix not factorised')
    call check(m == b%m,'LU_backsub_complex_mat','rhs has the wrong dimensions')
    call zgetrs	(trans_char,m,nrhs,self%elem,m,self%ipiv,b%elem,m,info)
    call check(info == 0,'LU_backsub_complex_mat','problem with back substitution')
    return
  end subroutine LU_backsub_complex_mat


  
  !---------------------------------------!
  !   constructors for complex matrices   !
  !---------------------------------------!

  type(complex_mat) function allocate_complex_mat_from_dimensions(m,n) result(a)
    integer(i4b), intent(in) :: m,n
    a%m = m
    a%n = n
    allocate(a%elem(m,n))
    a%elem = 0.0_dp
    if(m == n) a%square = .true.
    a%allocated = .true.
    return
  end function allocate_complex_mat_from_dimensions


  type(complex_mat) function allocate_square_complex_mat_from_dimensions(m) result(a)
    integer(i4b), intent(in) :: m    
    a = complex_mat(m,m)
    return
  end function allocate_square_complex_mat_from_dimensions


  type(complex_mat) function allocate_complex_mat_from_complex_mat(b) result(a)
    type(complex_mat), intent(in) :: b
    a%m = b%m
    a%n = b%n
    allocate(a%elem,source = b%elem)
    a%square = b%square
    a%allocated = .true.
    return
  end function allocate_complex_mat_from_complex_mat


  type(complex_mat) function allocate_complex_mat_from_rmat(b) result(a)
    type(rmat), intent(in) :: b
    a%m = b%m
    a%n = b%n
    allocate(a%elem(a%m,a%n))
    a%elem = b%elem
    a%square = b%square
    a%allocated = .true.
    return
  end function allocate_complex_mat_from_rmat


  type(complex_mat) function allocate_complex_mat_from_array(b) result(a)
    complex(dpc), dimension(:,:), intent(in) :: b
    a%m = size(b,1)
    a%n = size(b,2)
    if(a%m == a%n) a%square = .true.
    allocate(a%elem,source = b)
    a%allocated = .true.
    return
  end function allocate_complex_mat_from_array


  type(complex_mat) function unit_complex_mat(m) result(a)
    integer(i4b), intent(in) :: m
    integer(i4b) :: i
    a = complex_mat(m)
    do i = 1,m
       call a%inc(i,i,1.0_dp)
    end do
    return
  end function unit_complex_mat



  

  
  
end module module_LAPACK
