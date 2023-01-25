module module_LAPACK

  use module_constants
  use module_error
  implicit none


  interface allocate_vector
     procedure :: allocate_vector_real
     procedure :: allocate_vector_complex     
  end interface allocate_vector

  interface allocate_matrix
     procedure :: allocate_matrix_real_main
     procedure :: allocate_matrix_real_reduced
     procedure :: allocate_matrix_complex_main
     procedure :: allocate_matrix_complex_reduced
  end interface allocate_matrix
  
  interface allocate_matrix_bs
     procedure :: allocate_matrix_bs_real
     procedure :: allocate_matrix_bs_complex
  end interface allocate_matrix_bs

  interface allocate_matrix_b
     procedure :: allocate_matrix_b_real_main
     procedure :: allocate_matrix_b_real_reduced
     procedure :: allocate_matrix_b_complex_main
     procedure :: allocate_matrix_b_complex_reduced
  end interface allocate_matrix_b

  interface row_index_matrix_b
     procedure :: row_index_matrix_b_main
     procedure :: row_index_matrix_b_reduced     
  end interface row_index_matrix_b

  interface cholesky_matrix_bs     
     procedure :: cholesky_matrix_bs_real
     procedure :: cholesky_matrix_bs_complex
  end interface cholesky_matrix_bs


  interface cholesky_backsub_matrix_bs
     procedure :: cholesky_backsub_matrix_bs_real_single_overwrite
     procedure :: cholesky_backsub_matrix_bs_real_single
     procedure :: cholesky_backsub_matrix_bs_real_many_overwrite
     procedure :: cholesky_backsub_matrix_bs_real_many
     procedure :: cholesky_backsub_matrix_bs_complex_single_overwrite
     procedure :: cholesky_backsub_matrix_bs_complex_single
     procedure :: cholesky_backsub_matrix_bs_complex_many_overwrite
     procedure :: cholesky_backsub_matrix_bs_complex_many     
  end interface cholesky_backsub_matrix_bs

  
contains

  !--------------------------------------------------------------------!
  !                         routines for vectors                       !
  !--------------------------------------------------------------------!

  subroutine allocate_vector_real(n,x)
    integer(i4b), intent(in) :: n
    real(dp), dimension(:), allocatable, intent(inout) :: x
    if(allocated(x)) then
       if(size(x) /= n) then
          deallocate(x)
          allocate(x(n))
       end if
    else
       allocate(x(n))
    end if
    x = 0.0_dp
    return
  end subroutine allocate_vector_real


  subroutine allocate_vector_complex(n,x)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(:), allocatable, intent(inout) :: x
    if(allocated(x)) then
       if(size(x) /= n) then
          deallocate(x)
          allocate(x(n))
       end if
    else
       allocate(x(n))
    end if
    x = 0.0_dp
    return
  end subroutine allocate_vector_complex



  !--------------------------------------------------------------------!
  !                         routines for matrices                       !
  !--------------------------------------------------------------------!
  
  subroutine allocate_matrix_real_main(m,n,a)
    integer(i4b), intent(in) :: m,n
    real(dp), dimension(:,:), allocatable, intent(inout) :: a
    if(allocated(a)) then
       if(size(a,1) /= m .or. size(a,2) /= n) then
          deallocate(a)
          allocate(a(m,n))
       end if
    else
       allocate(a(m,n))
    end if
    a = 0.0_dp
    return
  end subroutine allocate_matrix_real_main


  subroutine allocate_matrix_real_reduced(n,a)
    integer(i4b), intent(in) :: n
    real(dp), dimension(:,:), allocatable, intent(inout) :: a
    call allocate_matrix(n,n,a)
    return
  end subroutine allocate_matrix_real_reduced


  subroutine allocate_matrix_complex_main(m,n,a)
    integer(i4b), intent(in) :: m,n
    complex(dpc), dimension(:,:), allocatable, intent(inout) :: a
    if(allocated(a)) then
       if(size(a,1) /= m .or. size(a,2) /= n) then
          deallocate(a)
          allocate(a(m,n))
       end if
    else
       allocate(a(m,n))
    end if
    a = 0.0_dp
    return
  end subroutine allocate_matrix_complex_main


  subroutine allocate_matrix_complex_reduced(n,a)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(:,:), allocatable, intent(inout) :: a
    call allocate_matrix(n,n,a)
    return
  end subroutine allocate_matrix_complex_reduced
  

  
  
  !--------------------------------------------------------------------!
  !               routines for banded symmetric matrices               !
  !--------------------------------------------------------------------!
  

  subroutine allocate_matrix_bs_real(n,kd,a)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(:,:), allocatable, intent(inout) :: a
    integer(i4b) :: ld
    ld = kd+1
    if(allocated(a)) then
       if(size(a,1) /= ld .or. size(a,2) /= n) then
          deallocate(a)
          allocate(a(ld,n))
       end if
    else
       allocate(a(ld,n))
    end if
    a = 0.0_dp       
    return
  end subroutine allocate_matrix_bs_real

  subroutine allocate_matrix_bs_complex(n,kd,a)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(:,:), allocatable, intent(inout) :: a
    integer(i4b) :: ld
    ld = kd+1
    if(allocated(a)) then
       if(size(a,1) /= ld .or. size(a,2) /= n) then
          deallocate(a)
          allocate(a(ld,n))
       end if
    else
       allocate(a(ld,n))
    end if
    a = 0.0_dp       
    return
  end subroutine allocate_matrix_bs_complex

  integer(i4b) pure function row_index_matrix_bs(kd,i,j) result(k)
    integer(i4b), intent(in) :: kd,i,j
    k = kd+1+i-j
    return
  end function row_index_matrix_bs


  subroutine cholesky_matrix_bs_real(kd,a)
    integer(i4b), intent(in) :: kd
    real(dp), dimension(:,:), intent(inout) :: a
    integer(i4b) :: n,ld,info
    ld = size(a,1)
    call check(ld == kd+1,'cholesky_matrix_bs_real','kd inconsistent with ld')
    n = size(a,2)
    call dpbtrf	('U',n,kd,a,ld,info)
    call check(info == 0,'cholesky_matrix_bs_real','problem with decompsition')
    return
  end subroutine cholesky_matrix_bs_real


  subroutine cholesky_matrix_bs_complex(kd,a)
    integer(i4b), intent(in) :: kd
    complex(dpc), dimension(:,:), intent(inout) :: a
    integer(i4b) :: n,ld,info
    ld = size(a,1)
    call check(ld == kd+1,'cholesky_matrix_bs_complex','kd inconsistent with ld')
    n = size(a,2)
    call zpbtrf	('U',n,kd,a,ld,info)
    call check(info == 0,'cholesky_matrix_bs_complex','problem with decompsition')
    return
  end subroutine cholesky_matrix_bs_complex

  
  subroutine cholesky_backsub_matrix_bs_real_single_overwrite(kd,a,b)
    integer(i4b), intent(in) :: kd
    real(dp), dimension(:,:), intent(inout) :: a
    real(dp), dimension(:), intent(inout), target :: b
    integer(i4b) :: n,ld,info
    real(dp), dimension(:,:), pointer :: bl
    ld = size(a,1)
    call check(ld == kd+1,'cholesky_backsub_matrix_bs_real_single_overwrite', &
                          'kd inconsistent with ld')
    n = size(a,2)
    call check(size(b) == n,'cholesky_backsub_matrix_bs_real_single_overwrite', &
                            'b has the wrong dimension')
    bl(1:n,1:1) => b
    call dpbtrs('U',n,kd,1,a,ld,bl,n,info)
    call check(info == 0,'cholesky_backsub_matrix_bs_real_single_overwrite', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine cholesky_backsub_matrix_bs_real_single_overwrite


  subroutine cholesky_backsub_matrix_bs_real_single(kd,a,b,c)
    integer(i4b), intent(in) :: kd
    real(dp), dimension(:,:), intent(inout) :: a
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(:), intent(out), target :: c
    c = b
    call cholesky_backsub_matrix_bs(kd,a,c)
    return
  end subroutine cholesky_backsub_matrix_bs_real_single


  subroutine cholesky_backsub_matrix_bs_real_many_overwrite(kd,a,b)
    integer(i4b), intent(in) :: kd
    real(dp), dimension(:,:), intent(inout) :: a
    real(dp), dimension(:,:), intent(inout) :: b
    integer(i4b) :: n,ld,info,nrhs
    ld = size(a,1)
    call check(ld == kd+1,'cholesky_matrix_bs_real','kd inconsistent with ld')
    n = size(a,2)    
    call check(size(b,1) == n,'cholesky_backsub_matrix_bs_real_single_overwrite', &
                              'b has the wrong dimension')
    nrhs = size(b,2)
    call dpbtrs('U',nrhs,kd,1,a,ld,b,n,info)
    call check(info == 0,'cholesky_backsub_matrix_bs_real_many_overwrite', &
                         'problem with backsubtitution')
    return
  end subroutine cholesky_backsub_matrix_bs_real_many_overwrite


  subroutine cholesky_backsub_matrix_bs_real_many(kd,a,b,c)
    integer(i4b), intent(in) :: kd
    real(dp), dimension(:,:), intent(inout) :: a
    real(dp), dimension(:,:), intent(in)    :: b
    real(dp), dimension(:,:), intent(out)   :: c
    integer(i4b) :: n,nrhs
    n = size(b,1)
    nrhs = size(b,2)
    call check(size(c,1) == n,'cholesky_backsub_matrix_bs_real_single_overwrite', &
                              'b and c need same row dimension')
    call check(size(c,2) == nrhs,'cholesky_backsub_matrix_bs_real_single_overwrite', &
                                 'b and c need same column dimension')
    c = b
    call cholesky_backsub_matrix_bs(kd,a,c)
    return
  end subroutine cholesky_backsub_matrix_bs_real_many


  subroutine cholesky_backsub_matrix_bs_complex_single_overwrite(kd,a,b)
    integer(i4b), intent(in) :: kd
    complex(dpc), dimension(:,:), intent(inout) :: a
    complex(dpc), dimension(:), intent(inout), target :: b
    integer(i4b) :: n,ld,info
    complex(dpc), dimension(:,:), pointer :: bl
    ld = size(a,1)
    call check(ld == kd+1,'cholesky_backsub_matrix_bs_complex_single_overwrite', &
                          'kd inconsistent with ld')
    n = size(a,2)
    call check(size(b) == n,'cholesky_backsub_matrix_bs_complex_single_overwrite', &
                            'b has the wrong dimension')
    bl(1:n,1:1) => b
    call zpbtrs('U',n,kd,1,a,ld,bl,n,info)
    call check(info == 0,'cholesky_backsub_matrix_bs_complex_single_overwrite', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine cholesky_backsub_matrix_bs_complex_single_overwrite


  subroutine cholesky_backsub_matrix_bs_complex_single(kd,a,b,c)
    integer(i4b), intent(in) :: kd
    complex(dpc), dimension(:,:), intent(inout) :: a
    complex(dpc), dimension(:), intent(in) :: b
    complex(dpc), dimension(:), intent(out), target :: c
    c = b
    call cholesky_backsub_matrix_bs(kd,a,c)
    return
  end subroutine cholesky_backsub_matrix_bs_complex_single


  subroutine cholesky_backsub_matrix_bs_complex_many_overwrite(kd,a,b)
    integer(i4b), intent(in) :: kd
    complex(dpc), dimension(:,:), intent(inout) :: a
    complex(dpc), dimension(:,:), intent(inout) :: b
    integer(i4b) :: n,ld,info,nrhs
    ld = size(a,1)
    call check(ld == kd+1,'cholesky_backsub_matrix_bs_complex_many_overwrite', &
                          'kd inconsistent with ld')
    n = size(a,2)    
    call check(size(b,1) == n,'cholesky_backsub_matrix_bs_complex_single_overwrite', &
                              'b has the wrong dimension')
    nrhs = size(b,2)
    call zpbtrs('U',nrhs,kd,1,a,ld,b,n,info)
    call check(info == 0,'cholesky_backsub_matrix_bs_complex_many_overwrite', &
                         'problem with backsubtitution')
    return
  end subroutine cholesky_backsub_matrix_bs_complex_many_overwrite


  subroutine cholesky_backsub_matrix_bs_complex_many(kd,a,b,c)
    integer(i4b), intent(in) :: kd
    complex(dpc), dimension(:,:), intent(inout) :: a
    complex(dpc), dimension(:,:), intent(in)    :: b
    complex(dpc), dimension(:,:), intent(out)   :: c
    integer(i4b) :: n,nrhs
    n = size(b,1)
    nrhs = size(b,2)
    call check(size(c,1) == n,'cholesky_backsub_matrix_bs_complex_single_overwrite', &
                              'b and c need same row dimension')
    call check(size(c,2) == nrhs,'cholesky_backsub_matrix_bs_complex_single_overwrite', &
                                 'b and c need same column dimension')
    c = b
    call cholesky_backsub_matrix_bs(kd,a,c)
    return
  end subroutine cholesky_backsub_matrix_bs_complex_many


  
  !--------------------------------------------------------------------!
  !                    routines for banded matrices                    !
  !--------------------------------------------------------------------!

  
  subroutine allocate_matrix_b_real_main(m,n,kl,ku,a)
    integer(i4b), intent(in) :: m,n,kl,ku
    real(dp), dimension(:,:), allocatable, intent(inout) :: a
    integer(i4b) :: ld
    ld = 2*kl+ku+1
    if(allocated(a)) then
       if(size(a,1) /= ld .or. size(a,2) /= n) then
          deallocate(a)
          allocate(a(ld,n))
       end if
    else
       allocate(a(ld,n))
    end if
    a = 0.0_dp       
    return
  end subroutine allocate_matrix_b_real_main


  subroutine allocate_matrix_b_real_reduced(n,kd,a)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(:,:), allocatable, intent(inout) :: a
    call allocate_matrix_b(n,n,kd,kd,a)
    return
  end subroutine allocate_matrix_b_real_reduced


  subroutine allocate_matrix_b_complex_main(m,n,kl,ku,a)
    integer(i4b), intent(in) :: m,n,kl,ku
    complex(dpc), dimension(:,:), allocatable, intent(inout) :: a
    integer(i4b) :: ld
    ld = 2*kl+ku+1
    if(allocated(a)) then
       if(size(a,1) /= ld .or. size(a,2) /= n) then
          deallocate(a)
          allocate(a(ld,n))
       end if
    else
       allocate(a(ld,n))
    end if
    a = 0.0_dp       
    return
  end subroutine allocate_matrix_b_complex_main


  subroutine allocate_matrix_b_complex_reduced(n,kd,a)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(:,:), allocatable, intent(inout) :: a
    call allocate_matrix_b(n,n,kd,kd,a)
    return
  end subroutine allocate_matrix_b_complex_reduced
  
 
  integer(i4b) pure function row_index_matrix_b_main(kl,ku,i,j) result(k)
    integer(i4b), intent(in) :: kl,ku,i,j
    k = kl+ku+1+i-j
    return
  end function row_index_matrix_b_main


  integer(i4b) pure function row_index_matrix_b_reduced(kd,i,j) result(k)
    integer(i4b), intent(in) :: kd,i,j
    k = 2*kd+1+i-j
    return
  end function row_index_matrix_b_reduced
  

  !------------------------------------------------------------------------!
  !             routines to convert between matrix formats                 !
  !------------------------------------------------------------------------!


  subroutine matrix_real_to_matrix_complex(a,b)
    real(dp), dimension(:,:), intent(in) :: a
    complex(dpc), dimension(:,:), allocatable, intent(inout) :: b
    integer(i4b) :: m,n
    m = size(a,1)
    n = size(a,2)
    call allocate_matrix(m,n,b)
    b = a
    return
  end subroutine matrix_real_to_matrix_complex

  subroutine matrix_real_bs_to_matrix_real(kd,a,b)
    integer(i4b), intent(in) :: kd
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(:,:), allocatable, intent(inout) :: b
    integer(i4b) :: ld,n,i,j,k
    ld = size(a,1)
    n  = size(a,2)
    call check(ld == kd+1,'matrix_real_bs_to_matrix_real','ld and kd are inconsistent')
    call allocate_matrix(n,n,b)
    do j = 1,n
       k = 0
       do i = max(1,j-kd),j
          k = k+1
          b(i,j) = a(k,j)
          b(j,k) = a(k,j)
       end do
    end do
    return
  end subroutine matrix_real_bs_to_matrix_real


  subroutine matrix_complex_bs_to_matrix_complex(kd,a,b)
    integer(i4b), intent(in) :: kd
    complex(dpc), dimension(:,:), intent(in) :: a
    complex(dpc), dimension(:,:), allocatable, intent(inout) :: b
    integer(i4b) :: ld,n,i,j,k
    ld = size(a,1)
    n  = size(a,2)
    call check(ld == kd+1,'matrix_complex_bs_to_matrix_complex','ld and kd are inconsistent')
    call allocate_matrix(n,n,b)
    do j = 1,n
       k = 0
       do i = max(1,j-kd),j
          k = k+1
          b(i,j) = a(k,j)
          b(j,k) = a(k,j)
       end do
    end do
    return
  end subroutine matrix_complex_bs_to_matrix_complex
  
end module module_LAPACK
