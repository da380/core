module module_LAPACK

  use module_constants
  use module_error
  implicit none

  type real_matrix
     logical :: allocated = .false.
     integer(i4b) :: m = 0
     integer(i4b) :: n = 0
     real(dp), dimension(:,:), allocatable :: elem
   contains
     procedure :: deallocate => deallocate_real_matrix
  end type real_matrix

  interface real_matrix
     procedure :: allocate_real_matrix
     procedure :: allocate_diagonal_real_matrix
     procedure :: allocate_move_real_matrix
     procedure :: allocate_from_array_real_matrix
  end interface real_matrix
  
  type, extends(real_matrix) :: LU_real_matrix
     logical :: factorised = .false.
     integer(i4b), dimension(:), allocatable :: ipiv
   contains
     procedure :: deallocate => deallocate_LU_real_matrix
     procedure :: factorise  => factorise_LU_real_matrix
     procedure :: backsub => backsub_LU_real_matrix
  end type LU_real_matrix

  interface LU_real_matrix
     procedure :: allocate_LU_real_matrix
     procedure :: allocate_diagonal_LU_real_matrix
     procedure :: allocate_from_real_matrix_LU_real_matrix
  end interface LU_real_matrix


  type, extends(real_matrix) :: banded_real_matrix
     integer(i4b) :: kl
     integer(i4b) :: ku
     integer(i4b) :: ldab
   contains
     procedure :: deallocate => deallocate_banded_real_matrix
  end type banded_real_matrix
  
  
contains


  !===============================================!
  !          procedures for real matrices         !
  !===============================================!
  
  subroutine deallocate_real_matrix(self)
    class(real_matrix), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%elem)
    self%m = 0
    self%n = 0
    self%allocated = .false.    
    return
  end subroutine deallocate_real_matrix

  type(real_matrix) function allocate_real_matrix(m,n) result(a)
    integer(i4b), intent(in) :: m,n
    a%m = m
    a%n = n
    allocate(a%elem(m,n))
    a%elem = 0.0_dp
    a%allocated = .true.
    return
  end function allocate_real_matrix


  type(real_matrix) function allocate_diagonal_real_matrix(m) result(a)
    integer(i4b), intent(in) :: m
    a = real_matrix(m,m)
    return
  end function allocate_diagonal_real_matrix


  type(real_matrix) function allocate_move_real_matrix(b,save) result(a)
    type(real_matrix), intent(inout) :: b
    logical, intent(in), optional :: save
    logical :: save_local
    if(present(save)) then
       save_local = save
    else
       save_local = .false.
    end if
    a%m = b%m
    a%n = b%n
    if(save_local) then
       allocate(a%elem,source = b%elem)
    else
       call move_alloc(b%elem,a%elem)
       b%allocated = .false.
    end if
    a%allocated = .true.
    return
  end function allocate_move_real_matrix


  type(real_matrix) function allocate_from_array_real_matrix(b) result(a)
    real(dp), dimension(:,:), intent(in) :: b
    a%m = size(b,1)
    a%n = size(b,2)
    allocate(a%elem,source = b)
    a%allocated = .true.
    return
  end function allocate_from_array_real_matrix

  
  
  

  subroutine deallocate_LU_real_matrix(self)
    class(LU_real_matrix), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%ipiv)
    self%factorised = .false.
    call self%real_matrix%deallocate()    
    return
  end subroutine deallocate_LU_real_matrix

  subroutine factorise_LU_real_matrix(self)
    class(LU_real_matrix), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'factorise_LU_real_matrix','matrix not allocated')
    if(self%factorised) return
    call dgetrf	(self%m,self%n,self%elem,self%m,self%ipiv,info)
    call check(info == 0,'factorise_LU_real_matrix','problem with factorisation')
    self%factorised = .true.
    return
  end subroutine factorise_LU_real_matrix


  subroutine backsub_LU_real_matrix(self,rhs,trans)
    class(LU_real_matrix), intent(in) :: self
    type(real_matrix), intent(inout) :: rhs
    logical, intent(in), optional :: trans
    character(len=1) :: trans_char
    integer(i4b) :: info
    trans_char = 'N'
    if(present(trans)) then
       if(trans) trans_char = 'T'
    end if          
    call check(self%m == self%n,'backsub_single_LU_real_matrix','matrix needs to be square')
    call check(self%factorised,'backsub_single_LU_real_matrix','matrix needs to be factorised')
    call check(self%n == rhs%m,'backsub_single_LU_real_matrix','rhs wrong dimension')    
    call dgetrs	(trans_char,self%n,rhs%n,self%elem,self%n,self%ipiv,rhs%elem,self%n,info)
    call check(info == 0,'backsub_single_LU_real_matrix','problem with backsubstitution')
    return
  end subroutine backsub_LU_real_matrix




  type(LU_real_matrix) function allocate_LU_real_matrix(m,n) result(a)
    integer(i4b) :: m,n
    a%real_matrix = real_matrix(m,n)
    allocate(a%ipiv(min(m,n)))
    return
  end function allocate_LU_real_matrix

  type(LU_real_matrix) function allocate_diagonal_LU_real_matrix(m) result(a)
    integer(i4b) :: m
    a = LU_real_matrix(m,m)
    return
  end function allocate_diagonal_LU_real_matrix  

  
  type(LU_real_matrix) function allocate_from_real_matrix_LU_real_matrix(b,save) result(a)
    type(real_matrix), intent(inout) :: b
    logical, intent(in), optional :: save
    logical :: save_local 
    a%real_matrix = real_matrix(b,save)
    allocate(a%ipiv(min(a%m,a%n)))
    call a%factorise()
    return
  end function allocate_from_real_matrix_LU_real_matrix


  type(LU_real_matrix) function allocate_from_array_LU_real_matrix(b) result(a)
    real(dp), dimension(:,:), intent(in) :: b
    a = LU_real_matrix(real_matrix(b))
    return
  end function allocate_from_array_LU_real_matrix


  !======================================================!
  !          procedures for banded real matrices         !
  !======================================================!
  
end module module_LAPACK
