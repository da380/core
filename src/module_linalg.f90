module module_linalg


  use module_constants
  use module_error
  implicit none

  type :: matrix_list

     logical :: allocated = .false.
     logical :: real = .true.     
     integer(i4b) :: m = 0
     integer(i4b) :: n = 0     
     integer(i4b) :: saved = 0
     integer(i4b) :: avail = 0
     integer(i4b), dimension(:), allocatable :: row
     integer(i4b), dimension(:), allocatable :: col     
     real(dp), dimension(:), allocatable :: rdat
     complex(dpc), dimension(:), allocatable :: cdat

   contains

     procedure :: deallocate => deallocate_matrix_list
     procedure :: reallocate => reallocate_matrix_list
     procedure :: print => print_matrix_list
     procedure, private :: add_real_to_matrix_list
     procedure, private :: add_real_block_to_matrix_list
     procedure, private :: add_complex_to_matrix_list
     generic, public :: add => add_real_to_matrix_list,       &
                               add_real_block_to_matrix_list, &
                               add_complex_to_matrix_list
  end type matrix_list

  interface matrix_list
     procedure :: allocate_matrix_list
  end interface matrix_list
  

contains

  type(matrix_list) function allocate_matrix_list(m,n,real,avail) result(a)
    integer(i4b), intent(in) :: m
    integer(i4b), intent(in) :: n
    logical, intent(in), optional :: real
    integer(i4b), intent(in), optional :: avail
    if(present(real)) a%real = real
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
    if(allocated(self%col)) deallocate(self%col)
    if(allocated(self%row)) deallocate(self%row)
    self%avail = 0
    self%saved = 0
    self%n = 0
    self%m = 0
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

    integer(i4b) :: k

    ! check if reallocatation is needed
    k = self%saved + 1
    if(k > self%avail) call self%reallocate()   
    
    self%row(k) = i
    self%col(k) = j
    if(self%real) then
       self%rdat(k) = val
    else
       self%cdat(k) = val
    end if
    self%saved = k
    
    return
  end subroutine add_real_to_matrix_list


  subroutine add_real_block_to_matrix_list(self,i,n,j,m,vals)
    class(matrix_list), intent(inout) :: self
    integer(i4b), intent(in) :: i,j,n,m
    real(dp), dimension(n,m), intent(in) :: vals
    
    integer(i4b) :: k,add,il,jl,ib,jb

    ! check if reallocation if needed
    k = k + m*n
    if(k > self%avail) then
       add = max(m*n,2*self%avail)
       call self%reallocate(add)
    end if

    k = self%saved
    do jb = 1,n
       jl = j+jb-1
       do ib = 1,m
          il = i+ib-1
          k = k + 1
          self%row(k) = il
          self%col(k) = jl
          if(self%real) then
             self%rdat(k) = vals(ib,jb)
          else
             self%cdat(k) = vals(ib,jb)             
          end if          
       end do
    end do
    self%saved = k
    
    return
  end subroutine add_real_block_to_matrix_list
  

  subroutine add_complex_to_matrix_list(self,i,j,val)
    class(matrix_list), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val

    integer(i4b) :: saved

    ! check if reallocation is needed
    saved = self%saved + 1
    if(saved > self%avail) call self%reallocate()
    
    ! store the new data
    self%row(saved) = i
    self%col(saved) = j
    if(self%real) then
       self%rdat(saved) = real(val,kind=dp)
    else
       self%cdat(saved) = val
    end if
    self%saved = saved
    
    return
  end subroutine add_complex_to_matrix_list



  

  
end module module_linalg
