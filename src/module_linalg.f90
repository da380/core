module module_linalg


  use module_constants
  use module_error
  implicit none

  type :: matrix_list
     
     private
     logical :: allocated = .false.
     logical :: real = .true.
     logical :: herm = .false.
     integer(i4b) :: m = 0
     integer(i4b) :: n = 0     
     integer(i4b) :: saved = 0
     integer(i4b) :: avail = 0
     integer(i4b), dimension(:), allocatable :: row
     integer(i4b), dimension(:), allocatable :: col     
     real(dp), dimension(:), allocatable :: rdat
     complex(dpc), dimension(:), allocatable :: cdat

   contains
     
     private
     procedure, public :: deallocate => deallocate_matrix_list
     procedure :: reallocate => reallocate_matrix_list
     procedure, public :: print => print_matrix_list
     procedure :: add_real_to_matrix_list
     procedure :: add_complex_to_matrix_list
     generic, public :: add => add_real_to_matrix_list,  & 
                        add_complex_to_matrix_list
     procedure :: matrix_list_to_real_dense
     procedure :: matrix_list_to_complex_dense
     generic, public :: dense => matrix_list_to_real_dense, &
                                 matrix_list_to_complex_dense
     procedure :: matrix_list_to_real_banded
     procedure :: matrix_list_to_complex_banded
     generic, public :: banded => matrix_list_to_real_banded, & 
                                  matrix_list_to_complex_banded
     procedure :: matrix_list_to_real_hermitian_banded
     procedure :: matrix_list_to_complex_hermitian_banded
     generic, public :: hermitian_banded => matrix_list_to_real_hermitian_banded, & 
                                  matrix_list_to_complex_hermitian_banded

     
     
  end type matrix_list

  interface matrix_list
     procedure :: allocate_matrix_list
  end interface matrix_list
  

contains

  type(matrix_list) function allocate_matrix_list(m,n,real,herm,avail) result(a)
    integer(i4b), intent(in) :: m
    integer(i4b), intent(in) :: n
    logical, intent(in), optional :: real
    logical, intent(in), optional :: herm
    integer(i4b), intent(in), optional :: avail
    a%m = m
    a%n = n
    if(present(real)) a%real = real
    if(present(herm)) then
       if(m == n) then
          a%herm = herm
       else
          print *, 'allocate_matrix_list: matrix dimensions do not match. Hermitian option ignored.'
       end if
    end if
    if(present(avail)) then
       a%avail = avail
    else
       a%avail = 5*max(m,n)
    end if
    allocate(a%row(a%avail))
    allocate(a%col(a%avail))    
    if(a%real) then
       allocate(a%rdat(a%avail))
    else
       allocate(a%cdat(a%avail))
    end if
    a%allocated = .true.
    return
  end function allocate_matrix_list


  subroutine deallocate_matrix_list(self)
    class(matrix_list), intent(inout) :: self
    if(allocated(self%cdat)) deallocate(self%cdat)
    if(allocated(self%rdat)) deallocate(self%rdat)
    if(allocated(self%col))  deallocate(self%col)
    if(allocated(self%row))  deallocate(self%row)
    self%avail = 0
    self%saved = 0
    self%n = 0
    self%m = 0
    self%herm = .false.
    self%real = .true.
    self%allocated = .false.
    return
  end subroutine deallocate_matrix_list

  
  subroutine reallocate_matrix_list(self,add)
    class(matrix_list), intent(inout) :: self
    integer(i4b), intent(in), optional :: add

    integer(i4b) :: avail_new,saved
    integer(i4b), dimension(:), allocatable :: row,col
    real(dp), dimension(:), allocatable :: rdat
    complex(dpc), dimension(:), allocatable :: cdat

    if(present(add)) then       
       avail_new = self%avail + add
       call check(avail_new >= 0,'reallocate_matrix_list','avail_new < 0')
    else
       avail_new = max(2*self%avail,1)
    end if
    saved = min(self%saved,avail_new)

    ! copy the row information
    allocate(row(avail_new))
    row(1:saved) = self%row(1:saved)
    call move_alloc(row,self%row)

    ! copy the column information
    allocate(col(avail_new))
    col(1:saved) = self%col(1:saved)
    call move_alloc(col,self%col)

    ! copy the data
    if(self%real) then
       allocate(rdat(avail_new))
       rdat(1:saved) = self%rdat(1:saved)
       call move_alloc(rdat,self%rdat)
    else
       allocate(cdat(avail_new))
       cdat(1:saved) = self%cdat(1:saved)
       call move_alloc(cdat,self%cdat)
    end if

    ! update the available space
    self%avail = avail_new    
    
    return
  end subroutine reallocate_matrix_list


  subroutine print_matrix_list(self)
    class(matrix_list), intent(in) :: self

    integer(i4b) :: k
    
    print *, ''
    print *, 'saved = ',self%saved,', available = ',self%avail
    print *, 'elements:'
    do k = 1,self%saved
       if(self%real) then
          print *, self%row(k),self%col(k),self%rdat(k)
       else
          print *, self%row(k),self%col(k),'(',real(self%cdat(k)),',',imag(self%cdat(k)),')'
       end if
    end do
    
    print *, ''
    
    return
  end subroutine print_matrix_list


  
  subroutine add_real_to_matrix_list(self,i,j,val)
    class(matrix_list), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val

    integer(i4b) :: k,il,jl

    ! check if reallocatation is needed
    k = self%saved + 1
    if(k > self%avail) call self%reallocate()   

    ! check for required symmetry
    if(self%herm .and. i > j) then
       il = j
       jl = i
    else
       il = i
       jl = j
    end if

    ! store the new data
    self%row(k) = il
    self%col(k) = jl
    if(self%real) then
       self%rdat(k) = val
    else
       self%cdat(k) = val
    end if
    self%saved = k
    
    return
  end subroutine add_real_to_matrix_list


  
  subroutine add_complex_to_matrix_list(self,i,j,val)
    class(matrix_list), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val

    integer(i4b) :: k,il,jl
    complex(dpc) :: lval

    ! check if reallocation is needed
    k= self%saved + 1
    if(k > self%avail) call self%reallocate()

    ! deal with required symmetry
    if(self%herm .and. i > j) then
       il = j
       jl = i
       lval = conjg(val)
    else
       il = i
       jl = j
       lval = val
    end if
    
    ! store the new data
    self%row(k) = il
    self%col(k) = jl
    if(self%real) then
       self%rdat(k) = real(lval,kind=dp)
    else
       self%cdat(k) = lval
    end if
    self%saved = k
    
    return
  end subroutine add_complex_to_matrix_list



  subroutine matrix_list_to_real_dense(self,a)
    class(matrix_list), intent(in) :: self
    real(dp), dimension(:,:), intent(out) :: a

    integer(i4b) :: i,j,k
    real(dp) :: val

    call check(size(a,1) == self%m,'matrix_list_to_real_dense','row dimensions do not match!')
    call check(size(a,2) == self%n,'matrix_list_to_real_dense','col dimensions do not match!')
        
    a = 0.0_dp
    do k = 1,self%saved
       i = self%row(k)
       j = self%col(k)
       if(self%real) then
          val = self%rdat(k)
       else
          val = real(self%cdat(k),kind=dp)
       end if
       a(i,j) = a(i,j) + val
       if(self%herm .and. i < j) then
          a(j,i) = a(j,i) + val
       end if
    end do
    
    return
  end subroutine matrix_list_to_real_dense



  subroutine matrix_list_to_complex_dense(self,a)
    class(matrix_list), intent(in) :: self
    complex(dpc), dimension(:,:), intent(out) :: a

    integer(i4b) :: i,j,k
    complex(dp) :: val

    call check(size(a,1) == self%m,'matrix_list_to_complex_dense','row dimensions do not match!')
    call check(size(a,2) == self%n,'matrix_list_to_complex_dense','col dimensions do not match!')
    
    a = 0.0_dp
    do k = 1,self%saved
       i = self%row(k)
       j = self%col(k)
       if(self%real) then
          val = self%rdat(k)
       else
          val = self%cdat(k)
       end if
       a(i,j) = a(i,j) + val
       if(self%herm .and. i < j) then
          a(j,i) = a(j,i) + conjg(val)
       end if
    end do
    
    return
  end subroutine matrix_list_to_complex_dense
  

  subroutine matrix_list_to_real_banded(self,kl,ku,a)
    class(matrix_list), intent(in) :: self
    integer(i4b), intent(in) :: kl,ku
    real(dp), dimension(:,:), intent(out) :: a

    integer(i4b) :: i,j,k,l,ld,imin,imax
    real(dp) :: val

    ! set the leading dimension
    ld = 2*kl+ku+1

    ! check the dimensions
    call check(size(a,1) == ld,    'matrix_list_to_real_banded','row dimension incompatible!')
    call check(size(a,2) == self%n,'matrix_list_to_real_banded','col dimensions do not match!')

    ! fill out the matrix
    a = 0.0_dp
    do k = 1,self%saved
       i = self%row(k)
       j = self%col(k)
       imin = max(1,j-ku)
       imax = min(self%n,k+kl)
       call check(i >= imin .and. i <= imax,'matrix_list_to_real_banded', &
                                            'element out of range')
       if(self%real) then
          val = self%rdat(k)
       else
          val = real(self%cdat(k),kind=dp)
       end if
       l = kl+ku+1+i-j
       a(l,j) = a(l,j) + val
       if(self%herm .and. i < j) then
          l = kl+ku+1+j-i
          a(l,i) = a(l,i) + val
       end if
    end do
    
    return
  end subroutine matrix_list_to_real_banded



  subroutine matrix_list_to_complex_banded(self,kl,ku,a)
    class(matrix_list), intent(in) :: self
    integer(i4b), intent(in) :: kl,ku
    complex(dpc), dimension(:,:), intent(out) :: a

    integer(i4b) :: i,j,k,l,ld,imin,imax
    complex(dpc) :: val

    ! set the leading dimension
    ld = 2*kl+ku+1

    ! check the dimensions
    call check(size(a,1) == ld,    'matrix_list_to_complex_banded','row dimension incompatible!')
    call check(size(a,2) == self%n,'matrix_list_to_complex_banded','col dimensions do not match!')

    ! fill out the matrix
    a = 0.0_dp
    do k = 1,self%saved
       i = self%row(k)
       j = self%col(k)
       imin = max(1,j-ku)
       imax = min(self%n,k+kl)
       call check(i >= imin .and. i <= imax,'matrix_list_to_complex_banded', &
                                            'element out of range')
       if(self%real) then
          val = self%rdat(k)
       else
          val = self%cdat(k)
       end if
       l = kl+ku+1+i-j
       a(l,j) = a(l,j) + val
       if(self%herm .and. i < j) then
          l = kl+ku+1+j-i
          a(l,i) = a(l,i) + conjg(val)
       end if
    end do
    
    return
  end subroutine matrix_list_to_complex_banded



  subroutine matrix_list_to_real_hermitian_banded(self,kd,a)
    class(matrix_list), intent(in) :: self
    integer(i4b), intent(in) :: kd
    real(dp), dimension(:,:), intent(out) :: a

    integer(i4b) :: i,j,k,l,ld,imin
    real(dp) :: val

    call check(self%herm,'matrix_list_to_real_hermitian_banded','matrix list not hermitian')
    
    ! set the leading dimension
    ld = kd+1

    ! check the dimensions
    call check(size(a,1) == ld,    'matrix_list_to_real_hermitian_banded','row dimension incompatible!')
    call check(size(a,2) == self%n,'matrix_list_to_real_hermitian_banded','col dimensions do not match!')

    ! fill out the matrix
    a = 0.0_dp
    do k = 1,self%saved
       i = self%row(k)
       j = self%col(k)
       imin = max(1,j-kd)
       call check(i >= imin .and. i <= j, 'matrix_list_to_real_hermitian_banded', &
                                          'element out of range')
       if(self%real) then
          val = self%rdat(k)
       else
          val = real(self%cdat(k),kind=dp)
       end if
       l = kd+1+i-j
       a(l,j) = a(l,j) + val
    end do
    
    return
  end subroutine matrix_list_to_real_hermitian_banded


  subroutine matrix_list_to_complex_hermitian_banded(self,kd,a)
    class(matrix_list), intent(in) :: self
    integer(i4b), intent(in) :: kd
    complex(dpc), dimension(:,:), intent(out) :: a

    integer(i4b) :: i,j,k,l,ld,imin
    complex(dpc) :: val

    call check(self%herm,'matrix_list_to_complex_hermitian_banded','matrix list not hermitian')
    
    ! set the leading dimension
    ld = kd+1

    ! check the dimensions
    call check(size(a,1) == ld,    'matrix_list_to_complex_hermitian_banded','row dimension incompatible!')
    call check(size(a,2) == self%n,'matrix_list_to_complex_hermitian_banded','col dimensions do not match!')

    ! fill out the matrix
    a = 0.0_dp
    do k = 1,self%saved
       i = self%row(k)
       j = self%col(k)
       imin = max(1,j-kd)
       call check(i >= imin .and. i <= j, 'matrix_list_to_complex_hermitian_banded', &
                                          'element out of range')
       if(self%real) then
          val = self%rdat(k)
       else
          val = self%cdat(k)
       end if
       l = kd+1+i-j
       a(l,j) = a(l,j) + val
    end do
    
    return
  end subroutine matrix_list_to_complex_hermitian_banded
  
  
end module module_linalg
