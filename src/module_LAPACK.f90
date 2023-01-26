module module_LAPACK

  use module_constants
  use module_error
  implicit none

  !------------------------------------------------------!
  !                    vector declarations               !
  !------------------------------------------------------!

  interface allocate_vector
     procedure :: allocate_vector_real
     procedure :: allocate_vector_complex     
  end interface allocate_vector

  !------------------------------------------------------!
  !                    matrix declarations               !
  !------------------------------------------------------!
  
  interface allocate_matrix
     procedure :: allocate_matrix_real
     procedure :: allocate_matrix_complex
  end interface allocate_matrix
  
  interface factorise_matrix
     procedure :: factorise_matrix_real
     procedure :: factorise_matrix_complex
  end interface factorise_matrix
  
  interface backsub_matrix
     procedure :: backsub_matrix_real_single
     procedure :: backsub_matrix_real_many
     procedure :: backsub_matrix_complex_single
     procedure :: backsub_matrix_complex_many
  end interface backsub_matrix


  !------------------------------------------------------!
  !             symmetric matrix declarations            !
  !------------------------------------------------------!
  
  interface allocate_matrix_s
     procedure :: allocate_matrix_s_real
     procedure :: allocate_matrix_s_complex
  end interface allocate_matrix_s
  
  interface factorise_matrix_s
     procedure :: factorise_matrix_s_real
     procedure :: factorise_matrix_s_complex
  end interface factorise_matrix_s
  
  interface backsub_matrix_s
     procedure :: backsub_matrix_s_real_single
     procedure :: backsub_matrix_s_real_many
     procedure :: backsub_matrix_s_complex_single
     procedure :: backsub_matrix_s_complex_many
  end interface backsub_matrix_s
  
  !------------------------------------------------------!
  !               banded matrix declarations             !
  !------------------------------------------------------!

  interface allocate_matrix_b
     procedure :: allocate_matrix_b_real
     procedure :: allocate_matrix_b_complex
  end interface allocate_matrix_b

  interface row_index_matrix_b
     procedure :: row_index_matrix_b
  end interface row_index_matrix_b

  interface factorise_matrix_b
     procedure :: factorise_matrix_b_real
     procedure :: factorise_matrix_b_complex
  end interface factorise_matrix_b
  
  interface backsub_matrix_b
     procedure :: backsub_matrix_b_real_single
     procedure :: backsub_matrix_b_real_many
     procedure :: backsub_matrix_b_complex_single
     procedure :: backsub_matrix_b_complex_many
  end interface backsub_matrix_b

  
  !------------------------------------------------------!
  !          banded symmetric matrix declarations        !
  !------------------------------------------------------!
  
  interface allocate_matrix_bs
     procedure :: allocate_matrix_bs_real
     procedure :: allocate_matrix_bs_complex
  end interface allocate_matrix_bs
  
  interface factorise_matrix_bs     
     procedure :: factorise_matrix_bs_real
     procedure :: factorise_matrix_bs_complex
  end interface factorise_matrix_bs

  interface backsub_matrix_bs
     procedure :: backsub_matrix_bs_real_single
     procedure :: backsub_matrix_bs_real_many
     procedure :: backsub_matrix_bs_complex_single
     procedure :: backsub_matrix_bs_complex_many
  end interface backsub_matrix_bs


  !-----------------------------------------------!
  !               other declarations              !
  !-----------------------------------------------!

  interface matrix_s_to_matrix
     procedure ::  matrix_s_to_matrix_real
     procedure ::  matrix_s_to_matrix_complex
  end interface matrix_s_to_matrix

  interface matrix_b_to_matrix
     procedure ::  matrix_b_to_matrix_real
     procedure ::  matrix_b_to_matrix_complex
  end interface matrix_b_to_matrix

  interface matrix_bs_to_matrix
     procedure ::  matrix_bs_to_matrix_real
     procedure ::  matrix_bs_to_matrix_complex
  end interface matrix_bs_to_matrix

  interface matrix_bs_to_matrix_s
     procedure ::  matrix_bs_to_matrix_s_real
     procedure ::  matrix_bs_to_matrix_s_complex
  end interface matrix_bs_to_matrix_s
  
  interface matrix_bs_to_matrix_b
     procedure ::  matrix_bs_to_matrix_b_real
     procedure ::  matrix_bs_to_matrix_b_complex
  end interface matrix_bs_to_matrix_b
  
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
  
  subroutine allocate_matrix_real(m,n,a)
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
  end subroutine allocate_matrix_real


  subroutine allocate_matrix_complex(m,n,a)
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
  end subroutine allocate_matrix_complex
  

  subroutine factorise_matrix_real(m,n,a,ipiv)
    integer(i4b), intent(in) :: m,n
    real(dp), dimension(m,n), intent(inout) :: a
    integer(i4b), dimension(:), intent(inout), allocatable :: ipiv
    integer(i4b) :: p,info
    p = min(m,n)
    if(allocated(ipiv)) then
       if(size(ipiv) /= p) then
          deallocate(ipiv)
          allocate(ipiv(p))
       end if
    else
       allocate(ipiv(p))
    end if
    call dgetrf(m,n,a,m,ipiv,info)
    call check(info == 0,'factorise_matrix_real','problem with decomposition')
    return
  end subroutine factorise_matrix_real


  subroutine factorise_matrix_complex(m,n,a,ipiv)
    integer(i4b), intent(in) :: m,n
    complex(dpc), dimension(m,n), intent(inout) :: a
    integer(i4b), dimension(:), intent(inout), allocatable :: ipiv
    integer(i4b) :: p,info
    p = min(m,n)
    if(allocated(ipiv)) then
       if(size(ipiv) /= p) then
          deallocate(ipiv)
          allocate(ipiv(p))
       end if
    else
       allocate(ipiv(p))
    end if
    call zgetrf(m,n,a,m,ipiv,info)
    call check(info == 0,'factorise_matrix_complex','problem with decomposition')
    return
  end subroutine factorise_matrix_complex

  
  subroutine backsub_matrix_real_single(n,a,ipiv,b,trans)
    integer(i4b), intent(in) :: n
    real(dp), dimension(n,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    real(dp), dimension(n), intent(inout), target :: b
    character(len=1), intent(in), optional  :: trans
    character(len=1) :: transl
    integer(i4b) :: info
    real(dp), dimension(:,:), pointer :: bl
    if(present(trans)) then
       transl = trans
    else
       transl = 'N'
    end if
    bl(1:n,1:1) => b
    call dgetrs(transl,n,1,a,n,ipiv,bl,n,info)
    call check(info == 0,'backsub_matrix_real_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_real_single


  subroutine backsub_matrix_real_many(n,a,ipiv,nrhs,b,trans)
    integer(i4b), intent(in) :: n
    real(dp), dimension(n,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    integer(i4b), intent(in) :: nrhs
    real(dp), dimension(n,nrhs), intent(inout) :: b
    character(len=1), intent(in), optional  :: trans
    character(len=1) :: transl
    integer(i4b) :: info
    if(present(trans)) then
       transl = trans
    else
       transl = 'N'
    end if
    call dgetrs(transl,n,nrhs,a,n,ipiv,b,n,info)
    call check(info == 0,'backsub_matrix_real_many', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_real_many


  subroutine backsub_matrix_complex_single(n,a,ipiv,b,trans)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    complex(dpc), dimension(n), intent(inout), target :: b
    character(len=1), intent(in), optional  :: trans
    character(len=1) :: transl
    integer(i4b) :: info
    complex(dpc), dimension(:,:), pointer :: bl
    if(present(trans)) then
       transl = trans
    else
       transl = 'N'
    end if
    bl(1:n,1:1) => b
    call zgetrs(transl,n,1,a,n,ipiv,bl,n,info)
    call check(info == 0,'backsub_matrix_complex_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_complex_single


  subroutine backsub_matrix_complex_many(n,a,ipiv,nrhs,b,trans)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    integer(i4b), intent(in) :: nrhs
    complex(dpc), dimension(n,nrhs), intent(inout) :: b
    character(len=1), intent(in), optional  :: trans
    character(len=1) :: transl
    integer(i4b) :: info
    if(present(trans)) then
       transl = trans
    else
       transl = 'N'
    end if
    call zgetrs(transl,n,nrhs,a,n,ipiv,b,n,info)
    call check(info == 0,'backsub_matrix_complex_many', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_complex_many


  !--------------------------------------------------------------------!
  !                   routines for symmetric matrices                  !
  !--------------------------------------------------------------------!

  subroutine allocate_matrix_s_real(n,a)
    integer(i4b), intent(in) :: n
    real(dp), dimension(:,:), allocatable :: a
    call allocate_matrix(n,n,a)
    return
  end subroutine allocate_matrix_s_real


  subroutine allocate_matrix_s_complex(n,a)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(:,:), allocatable :: a
    call allocate_matrix(n,n,a)
    return
  end subroutine allocate_matrix_s_complex


  subroutine factorise_matrix_s_real(n,a,ipiv)
    integer(i4b), intent(in) :: n
    real(dp), dimension(n,n) , intent(inout) :: a
    integer(i4b), dimension(:), allocatable, intent(inout) :: ipiv
    integer(i4b) :: lwork,info
    real(dp), dimension(1) :: work_tmp
    real(dp), dimension(:), allocatable :: work
    if(allocated(ipiv)) then
       if(size(ipiv) /= n) then
          deallocate(ipiv)
          allocate(ipiv(n))
       end if
    else
       allocate(ipiv(n))
    end if
    call dsytrf('U',n,a,n,ipiv,work_tmp,-1,info)
    lwork = work_tmp(1)
    allocate(work(lwork))
    call dsytrf('U',n,a,n,ipiv,work,lwork,info)
    call check(info == 0,'factorise_matrix_s_real','problem with factorisation')
    return
  end subroutine factorise_matrix_s_real


  subroutine factorise_matrix_s_complex(n,a,ipiv)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n,n) , intent(inout) :: a
    integer(i4b), dimension(:), allocatable, intent(inout) :: ipiv
    integer(i4b) :: lwork,info
    complex(dpc), dimension(1) :: work_tmp
    complex(dpc), dimension(:), allocatable :: work
    if(allocated(ipiv)) then
       if(size(ipiv) /= n) then
          deallocate(ipiv)
          allocate(ipiv(n))
       end if
    else
       allocate(ipiv(n))
    end if
    call zsytrf('U',n,a,n,ipiv,work_tmp,-1,info)
    lwork = work_tmp(1)
    allocate(work(lwork))
    call zsytrf('U',n,a,n,ipiv,work,lwork,info)
    call check(info == 0,'factorise_matrix_s_complex','problem with factorisation')
    return
  end subroutine factorise_matrix_s_complex


  subroutine backsub_matrix_s_real_single(n,a,ipiv,b)
    integer(i4b), intent(in) :: n
    real(dp), dimension(n,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    real(dp), dimension(n), intent(inout), target :: b
    integer(i4b) :: info
    real(dp), dimension(:,:), pointer :: bl
    bl(1:n,1:1) => b
    call dsytrs('U',n,1,a,n,ipiv,bl,n,info) 
    call check(info == 0,'backsub_matrix_real_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_s_real_single

 
  subroutine backsub_matrix_s_real_many(n,a,ipiv,nrhs,b)
    integer(i4b), intent(in) :: n
    real(dp), dimension(n,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    integer(i4b), intent(in) :: nrhs
    real(dp), dimension(n,nrhs), intent(inout) :: b
    integer(i4b) :: info
    call dsytrs('U',n,nrhs,a,n,ipiv,b,n,info) 
    call check(info == 0,'backsub_matrix_real_single', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_s_real_many


  subroutine backsub_matrix_s_complex_single(n,a,ipiv,b)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    complex(dpc), dimension(n), intent(inout), target :: b
    integer(i4b) :: info
    complex(dpc), dimension(:,:), pointer :: bl
    bl(1:n,1:1) => b
    call zsytrs('U',n,1,a,n,ipiv,bl,n,info) 
    call check(info == 0,'backsub_matrix_complex_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_s_complex_single

 
  subroutine backsub_matrix_s_complex_many(n,a,ipiv,nrhs,b)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    integer(i4b), intent(in) :: nrhs
    complex(dpc), dimension(n,nrhs), intent(inout) :: b
    integer(i4b) :: info
    call zsytrs('U',n,nrhs,a,n,ipiv,b,n,info) 
    call check(info == 0,'backsub_matrix_complex_single', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_s_complex_many


  !--------------------------------------------------------------------!
  !                    routines for banded matrices                    !
  !--------------------------------------------------------------------!

  
  subroutine allocate_matrix_b_real(m,n,kl,ku,a)
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
  end subroutine allocate_matrix_b_real


  subroutine allocate_matrix_b_complex(m,n,kl,ku,a)
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
  end subroutine allocate_matrix_b_complex
  
 
  integer(i4b) pure function row_index_matrix_b(kl,ku,i,j) result(k)
    integer(i4b), intent(in) :: kl,ku,i,j
    k = kl+ku+1+i-j
    return
  end function row_index_matrix_b
  

  subroutine factorise_matrix_b_real(m,n,kl,ku,a,ipiv)
    integer(i4b), intent(in) :: m,n,kl,ku
    real(dp), dimension(2*kl+ku+1,n), intent(inout) :: a
    integer(i4b), dimension(:), allocatable, intent(inout) :: ipiv
    integer(i4b) :: p,ld,info
    p = min(m,n)
    if(allocated(ipiv)) then
       if(size(ipiv) /= p) then
          deallocate(ipiv)
          allocate(ipiv(p))
       end if
    else
       allocate(ipiv(p))
    end if
    ld = 2*kl+ku+1
    call dgbtrf(m,n,kl,ku,a,ld,ipiv,info)
    call check(info == 0,'factorise_matrix_b_real','problem with decomposition')
    return
  end subroutine factorise_matrix_b_real


  subroutine factorise_matrix_b_complex(m,n,kl,ku,a,ipiv)
    integer(i4b), intent(in) :: m,n,kl,ku
    complex(dpc), dimension(2*kl+ku+1,n), intent(inout) :: a
    integer(i4b), dimension(:), allocatable, intent(inout) :: ipiv
    integer(i4b) :: ld,p,info
    p = min(m,n)
    if(allocated(ipiv)) then
       if(size(ipiv) /= p) then
          deallocate(ipiv)
          allocate(ipiv(p))
       end if
    else
       allocate(ipiv(p))
    end if
    ld = 2*kl+ku+1
    call zgbtrf(m,n,kl,ku,a,ld,ipiv,info)
    call check(info == 0,'factorise_matrix_b_complex','problem with decomposition')
    return
  end subroutine factorise_matrix_b_complex


  subroutine backsub_matrix_b_real_single(n,kl,ku,a,ipiv,b,trans)
    integer(i4b), intent(in) :: n,kl,ku
    real(dp), dimension(2*kl+ku+1,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    real(dp), dimension(n), intent(inout), target :: b
    character(len=1), intent(in), optional  :: trans
    character(len=1) :: transl
    integer(i4b) :: lda,info
    real(dp), dimension(:,:), pointer :: bl
    if(present(trans)) then
       transl = trans
    else
       transl = 'N'
    end if
    lda = 2*kl+ku+1
    bl(1:n,1:1) => b
    call dgbtrs(transl,n,kl,ku,1,a,lda,ipiv,bl,n,info)
    call check(info == 0,'backsub_matrix_b_real_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_b_real_single

  
  subroutine backsub_matrix_b_real_many(n,kl,ku,a,ipiv,nrhs,b,trans)
    integer(i4b), intent(in) :: n,kl,ku
    real(dp), dimension(2*kl+ku+1,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    integer(i4b), intent(in) :: nrhs
    real(dp), dimension(n,nrhs), intent(inout) :: b
    character(len=1), intent(in), optional  :: trans
    character(len=1) :: transl
    integer(i4b) :: lda,info
    if(present(trans)) then
       transl = trans
    else
       transl = 'N'
    end if
    lda = 2*kl+ku+1
    call dgbtrs(transl,n,kl,ku,nrhs,a,lda,ipiv,b,n,info)
    call check(info == 0,'backsub_matrix_b_real_many', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_b_real_many


  subroutine backsub_matrix_b_complex_single(n,kl,ku,a,ipiv,b,trans)
    integer(i4b), intent(in) :: n,kl,ku
    complex(dpc), dimension(2*kl+ku+1,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    complex(dpc), dimension(n), intent(inout), target :: b
    character(len=1), intent(in), optional  :: trans
    character(len=1) :: transl
    integer(i4b) :: lda,info
    complex(dpc), dimension(:,:), pointer :: bl
    if(present(trans)) then
       transl = trans
    else
       transl = 'N'
    end if
    lda = 2*kl+ku+1
    bl(1:n,1:1) => b
    call zgbtrs(transl,n,kl,ku,1,a,lda,ipiv,bl,n,info)
    call check(info == 0,'backsub_matrix_b_complex_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_b_complex_single

  
  subroutine backsub_matrix_b_complex_many(n,kl,ku,a,ipiv,nrhs,b,trans)
    integer(i4b), intent(in) :: n,kl,ku
    complex(dpc), dimension(2*kl+ku+1,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    integer(i4b), intent(in) :: nrhs
    complex(dpc), dimension(n,nrhs), intent(inout) :: b
    character(len=1), intent(in), optional  :: trans
    character(len=1) :: transl
    integer(i4b) :: lda,info
    if(present(trans)) then
       transl = trans
    else
       transl = 'N'
    end if
    lda = 2*kl+ku+1
    call zgbtrs(transl,n,kl,ku,nrhs,a,lda,ipiv,b,n,info)
    call check(info == 0,'backsub_matrix_b_complex_many', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_b_complex_many

  
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


  subroutine factorise_matrix_bs_real(n,kd,a)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(kd+1,n), intent(inout) :: a
    integer(i4b) :: ld,info
    ld = kd+1
    call dpbtrf	('U',n,kd,a,ld,info)
    call check(info == 0,'factorise_matrix_bs_real','problem with decompsition')
    return
  end subroutine factorise_matrix_bs_real


  subroutine factorise_matrix_bs_complex(n,kd,a)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(kd+1,n), intent(inout) :: a
    integer(i4b) :: ld,info
    ld = kd+1
    call zpbtrf	('U',n,kd,a,ld,info)
    call check(info == 0,'factorise_matrix_bs_complex','problem with decompsition')
    return
  end subroutine factorise_matrix_bs_complex

  
  subroutine backsub_matrix_bs_real_single(n,kd,a,b)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(kd+1,n), intent(inout) :: a
    real(dp), dimension(n), intent(inout), target :: b
    integer(i4b) :: ld,info
    real(dp), dimension(:,:), pointer :: bl
    ld = kd+1
    bl(1:n,1:1) => b
    call dpbtrs('U',n,kd,1,a,ld,bl,n,info)
    call check(info == 0,'backsub_matrix_bs_real_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_bs_real_single


  subroutine backsub_matrix_bs_real_many(n,kd,a,nrhs,b)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(kd+1,n), intent(inout) :: a
    integer(i4b), intent(in) :: nrhs
    real(dp), dimension(n,nrhs), intent(inout) :: b
    integer(i4b) :: ld,info
    ld = kd+1
    call dpbtrs('U',n,kd,nrhs,a,ld,b,n,info)
    call check(info == 0,'backsub_matrix_bs_real_many', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_bs_real_many


  subroutine backsub_matrix_bs_complex_single(n,kd,a,b)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(kd+1,n), intent(inout) :: a
    complex(dpc), dimension(n), intent(inout), target :: b
    integer(i4b) :: ld,info
    complex(dpc), dimension(:,:), pointer :: bl
    ld = kd+1
    bl(1:n,1:1) => b
    call zpbtrs('U',n,kd,1,a,ld,bl,n,info)
    call check(info == 0,'backsub_matrix_bs_complex_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_bs_complex_single


  subroutine backsub_matrix_bs_complex_many(n,kd,a,nrhs,b)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(kd+1,n), intent(inout) :: a
    integer(i4b), intent(in) :: nrhs
    complex(dpc), dimension(n,nrhs), intent(inout) :: b
    integer(i4b) :: ld,info
    ld = kd+1
    call zpbtrs('U',n,kd,nrhs,a,ld,b,n,info)
    call check(info == 0,'backsub_matrix_bs_complex_many', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_bs_complex_many


  !--------------------------------------------------------------------!
  !              routines to convert between storage types             !
  !--------------------------------------------------------------------!


  subroutine matrix_s_to_matrix_real(n,a,b)
    integer(i4b), intent(in) :: n
    real(dp), dimension(n,n), intent(in) :: a
    real(dp), dimension(n,n), intent(out) :: b
    integer(i4b) :: i,j
    do j = 1,n
       do i = 1,j-1
          b(i,j) = a(i,j)
          b(j,i) = a(i,j)
       end do
       b(j,j) = a(j,j)          
       end do
    return
  end subroutine matrix_s_to_matrix_real


  subroutine matrix_s_to_matrix_complex(n,a,b)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n,n), intent(in) :: a
    complex(dpc), dimension(n,n), intent(out) :: b    
    integer(i4b) :: i,j
    do j = 1,n
       do i = 1,j-1
          b(i,j) = a(i,j)
          b(j,i) = a(i,j)
       end do
       b(j,j) = a(j,j)          
       end do
    return
  end subroutine matrix_s_to_matrix_complex


  subroutine matrix_b_to_matrix_real(m,n,kl,ku,a,b)
    integer(i4b), intent(in) :: m,n,kl,ku
    real(dp), dimension(2*kl+ku+1,n), intent(in) :: a
    real(dp), dimension(m,n), intent(out) :: b
    integer(i4b) :: i,j,k
    b = 0.0_dp
    do j = 1,n
       do i= max(1,j-ku),min(m,j+kl)
          k = kl+ku+1+i-j
          b(i,j) = a(k,j)
       end do
    end do
    return
  end subroutine matrix_b_to_matrix_real


  subroutine matrix_b_to_matrix_complex(m,n,kl,ku,a,b)
    integer(i4b), intent(in) :: m,n,kl,ku
    complex(dpc), dimension(2*kl+ku+1,n), intent(in) :: a
    complex(dpc), dimension(m,n), intent(out) :: b
    integer(i4b) :: i,j,k
    b = 0.0_dp
    do j = 1,n
       do i= max(1,j-ku),min(m,j+kl)
          k = kl+ku+1+i-j
          b(i,j) = a(k,j)
       end do
    end do
    return
  end subroutine matrix_b_to_matrix_complex


  subroutine matrix_bs_to_matrix_real(n,kd,a,b)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(kd+1,n), intent(in) :: a
    real(dp), dimension(n,n), intent(out) :: b
    integer(i4b) :: i,j,k
    b = 0.0_dp
    do j =1,n
       do i = 1,max(1,j-kd),j-1
          k = kd+1+i-j
          b(i,j) = a(k,j)
          b(j,i) = a(k,j)
       end do
       k = kd+1
       b(j,j) = a(k,j)
    end do
    return
  end subroutine matrix_bs_to_matrix_real


  subroutine matrix_bs_to_matrix_complex(n,kd,a,b)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(kd+1,n), intent(in) :: a
    complex(dpc), dimension(n,n), intent(out) :: b
    integer(i4b) :: i,j,k
    b = 0.0_dp
    do j =1,n
       do i = 1,max(1,j-kd),j-1
          k = kd+1+i-j
          b(i,j) = a(k,j)
          b(j,i) = a(k,j)
       end do
       k = kd+1
       b(j,j) = a(k,j)
    end do
    return
  end subroutine matrix_bs_to_matrix_complex


  subroutine matrix_bs_to_matrix_s_real(n,kd,a,b)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(kd+1,n), intent(in) :: a
    real(dp), dimension(n,n), intent(out) :: b
    integer(i4b) :: i,j,k
    b = 0.0_dp
    do j =1,n
       do i = 1,max(1,j-kd),j
          k = kd+1+i-j
          b(i,j) = a(k,j) 
       end do
    end do
    return
  end subroutine matrix_bs_to_matrix_s_real


  subroutine matrix_bs_to_matrix_s_complex(n,kd,a,b)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(kd+1,n), intent(in) :: a
    complex(dpc), dimension(n,n), intent(out) :: b
    integer(i4b) :: i,j,k
    b = 0.0_dp
    do j =1,n
       do i = 1,max(1,j-kd),j
          k = kd+1+i-j
          b(i,j) = a(k,j)
       end do
    end do
    return
  end subroutine matrix_bs_to_matrix_s_complex
  

  subroutine matrix_bs_to_matrix_b_real(n,kd,a,kl,ku,b)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(kd+1,n), intent(in) :: a
    integer(i4b), intent(in) :: kl,ku
    real(dp), dimension(2*kl+ku+1,n), intent(out) :: b
    integer(i4b) :: i,j,k1,k2    
    call check(kl >= kd,'matrix_bs_to_matrix_b_real', &
                        'lower band width too small')
    call check(ku >= kd,'matrix_bs_to_matrix_b_real', &
                        'upper band width too small')
    b = 0.0_dp
    do j = 1,n
       do i = max(1,j-kd),j-1
          k1 = kd+1+i-j
          k2 = kl+ku+1+i-j
          b(k2,j) = a(k1,j)
          k2 = kl+ku+1+j-i
          b(k2,i) = a(k1,j)          
       end do
       k1 = kd+1
       k2 = kl+ku+1
       b(k2,j) = a(k1,j)
    end do
    return
  end subroutine matrix_bs_to_matrix_b_real


  subroutine matrix_bs_to_matrix_b_complex(n,kd,a,kl,ku,b)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(kd+1,n), intent(in) :: a
    integer(i4b), intent(in) :: kl,ku
    complex(dpc), dimension(2*kl+ku+1,n), intent(out) :: b
    integer(i4b) :: i,j,k1,k2    
    call check(kl >= kd,'matrix_bs_to_matrix_b_complex', &
                        'lower band width too small')
    call check(ku >= kd,'matrix_bs_to_matrix_b_complex', &
                        'upper band width too small')
    b = 0.0_dp
    do j = 1,n
       do i = max(1,j-kd),j-1
          k1 = kd+1+i-j
          k2 = kl+ku+1+i-j
          b(k2,j) = a(k1,j)
          k2 = kl+ku+1+j-i
          b(k2,i) = a(k1,j)          
       end do
       k1 = kd+1
       k2 = kl+ku+1
       b(k2,j) = a(k1,j)
    end do
    return
  end subroutine matrix_bs_to_matrix_b_complex




  
end module module_LAPACK
