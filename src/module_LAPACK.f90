 module module_LAPACK

  use module_constants
  use module_error
  implicit none
  
  type mat
     logical :: real = .true.
     logical :: allocated = .false.
     logical :: square = .false.
     logical :: check = .false.
     character(len=:), allocatable, private :: factor  
     integer(i4b) :: m = 0
     integer(i4b) :: n = 0
     real(dp), dimension(:,:), allocatable :: elem_r
     complex(dpc), dimension(:,:), allocatable :: elem_c
     integer(i4b), dimension(:), allocatable, private :: ipiv
   contains
     procedure :: deallocate => deallocate_mat
     procedure :: allocate => allocate_mat
     procedure :: set_r_mat
     procedure :: set_c_mat
     generic   :: set        => set_r_mat,set_c_mat
     procedure :: get_r_mat
     procedure :: get_c_mat
     generic   :: get        => get_r_mat,get_c_mat
     procedure :: inc_r_mat
     procedure :: inc_c_mat
     generic   :: inc        => inc_r_mat,inc_c_mat
     procedure :: set_r_sym_mat
     procedure :: set_c_sym_mat
     generic   :: set_sym    => set_r_sym_mat,set_c_sym_mat
     procedure :: inc_r_sym_mat
     procedure :: inc_c_sym_mat
     generic   :: inc_sym    => inc_r_sym_mat,inc_c_sym_mat
     procedure :: LU         => LU_mat
     procedure :: LU_bsub => LU_bsub_mat
     procedure :: fac        => factor_mat
     procedure :: bsub       => bsub_mat
  end type mat

  type, extends(mat) :: bmat
     logical :: band_set = .false.
     integer(i4b) :: kl
     integer(i4b) :: ku
     integer(i4b) :: ld
   contains
     procedure :: deallocate => deallocate_bmat
     procedure :: band => band_bmat
     procedure :: allocate => allocate_bmat
     procedure :: set_r_mat => set_r_bmat
     procedure :: set_c_mat => set_c_bmat
     procedure :: get_r_mat => get_r_bmat
     procedure :: get_c_mat => get_c_bmat
     procedure :: inc_r_mat => inc_r_bmat
     procedure :: inc_c_mat => inc_c_bmat
     procedure :: LU         => LU_bmat
     procedure :: LU_bsub => LU_bsub_bmat
     procedure :: row_ind    => row_ind_bmat
  end type bmat

  type, extends(mat) :: hbmat
     logical :: band_set = .false.
     integer(i4b) :: kd
     integer(i4b) :: ld
   contains
     procedure :: deallocate     => deallocate_hbmat
     procedure :: band => band_hbmat
     procedure :: allocate => allocate_hbmat
     procedure :: set_r_mat => set_r_hbmat
     procedure :: set_c_mat => set_c_hbmat
     procedure :: get_r_mat => get_r_hbmat
     procedure :: get_c_mat => get_c_hbmat
     procedure :: inc_r_mat => inc_r_hbmat
     procedure :: inc_c_mat => inc_c_hbmat
     procedure :: set_r_sym_mat => set_r_hbmat
     procedure :: set_c_sym_mat => set_c_hbmat
     procedure :: inc_r_sym_mat => inc_r_hbmat
     procedure :: inc_c_sym_mat => inc_c_hbmat
     procedure :: row_ind        => row_ind_hbmat
     procedure :: cholesky       => cholesky_hbmat
     procedure :: cholesky_bsub  => cholesky_bsub_hbmat
     procedure :: fac            => factor_hbmat
     procedure :: bsub           => bsub_hbmat
  end type hbmat

  
  
contains


  !===============================================!
  !          procedures for real matrices         !
  !===============================================!
  
  subroutine deallocate_mat(self)
    class(mat), intent(inout) :: self
    if(.not.self%allocated) return
    self%m = 0
    self%n = 0
    self%real = .true.
    if(allocated(self%ipiv)) deallocate(self%ipiv)
    if(allocated(self%elem_r)) deallocate(self%elem_r)
    if(allocated(self%elem_c)) deallocate(self%elem_c)
    if(allocated(self%factor)) deallocate(self%factor)
    self%allocated = .false.    
    return
  end subroutine deallocate_mat


  subroutine allocate_mat(self,m,n,real)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: m,n
    logical, intent(in), optional :: real
    if(present(real)) self%real = real
    self%m = m
    self%n = n
    if(self%real) then
       allocate(self%elem_r(m,n))
       self%elem_r = 0.0_dp
    else
       allocate(self%elem_c(m,n))
       self%elem_c = 0.0_dp
    end if
    if(m == n) self%square = .true.
    self%allocated = .true.
    return
  end subroutine allocate_mat

  
  subroutine set_r_mat(self,i,j,val)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'set_mat','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'set_mat','row index out of range')
       call check(j >= 1 .and. j <= self%n,'set_mat','column index out of range')
    end if
    if(self%real) then
       self%elem_r(i,j) = val
    else
       self%elem_c(i,j) = val
    end if
    return
  end subroutine set_r_mat


  subroutine set_c_mat(self,i,j,val)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'set_mat','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'set_mat','row index out of range')
       call check(j >= 1 .and. j <= self%n,'set_mat','column index out of range')
    end if
    if(self%real) then
       self%elem_r(i,j) = val
    else
       self%elem_c(i,j) = val
    end if
    return
  end subroutine set_c_mat


  subroutine get_r_mat(self,i,j,val)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(out) :: val
    if(self%check) then
       call check(self%allocated,'get_mat','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'get_mat','row index out of range')
       call check(j >= 1 .and. j <= self%n,'get_mat','column index out of range')
    end if
    if(self%real) then
       val  = self%elem_r(i,j)
    else
       val = self%elem_c(i,j)
    end if
    return
  end subroutine get_r_mat


  subroutine get_c_mat(self,i,j,val)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(out) :: val
    if(self%check) then
       call check(self%allocated,'get_mat','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'get_mat','row index out of range')
       call check(j >= 1 .and. j <= self%n,'get_mat','column index out of range')
    end if
    if(self%real) then
       val = self%elem_r(i,j) 
    else
       val = self%elem_c(i,j) 
    end if
    return
  end subroutine get_c_mat
  
  
  
  subroutine inc_r_mat(self,i,j,val)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'inc_mat','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'inc_mat','row index out of range')
       call check(j >= 1 .and. j <= self%n,'inc_mat','column index out of range')
    end if
    if(self%real) then
       self%elem_r(i,j) = self%elem_r(i,j) + val
    else
       self%elem_c(i,j) = self%elem_c(i,j) + val
    end if
    return
  end subroutine inc_r_mat


  subroutine inc_c_mat(self,i,j,val)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'inc_mat','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'inc_mat','row index out of range')
       call check(j >= 1 .and. j <= self%n,'inc_mat','column index out of range')
    end if
    if(self%real) then
       self%elem_r(i,j) = self%elem_r(i,j) + val
    else
       self%elem_c(i,j) = self%elem_c(i,j) + val
    end if
    return
  end subroutine inc_c_mat


  subroutine set_r_sym_mat(self,i,j,val)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    call self%set(i,j,val)
    if(i /= j) then
       call self%set(j,i,val)
    end if
    return
  end subroutine set_r_sym_mat


  subroutine set_c_sym_mat(self,i,j,val)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val
    call self%set(i,j,val)
    if(i /= j) then
       call self%set(j,i,val)
    end if
    return
  end subroutine set_c_sym_mat
  
  
  subroutine inc_r_sym_mat(self,i,j,val)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    call self%inc(i,j,val)
    if(i /= j) then
       call self%inc(j,i,val)
    end if
    return
  end subroutine inc_r_sym_mat


  subroutine inc_c_sym_mat(self,i,j,val)
    class(mat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val
    call self%inc(i,j,val)
    if(i /= j) then
       call self%inc(j,i,val)
    end if
    return
  end subroutine inc_c_sym_mat
  
    
  subroutine LU_mat(self)
    class(mat), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'LU_mat','matrix not allocated')
    if(self%factor == 'LU') return
    call check(self%factor == '','LU_mat','matrix already factored using different scheme')    
    allocate(self%ipiv(min(self%m,self%n)))
    if(self%real) then
       call dgetrf(self%m,self%n,self%elem_r,self%m,self%ipiv,info)
    else
       call zgetrf(self%m,self%n,self%elem_c,self%m,self%ipiv,info)
    end if
    call check(info == 0,'LU_mat','problem with factorisation')
    self%factor = 'LU'
    return
  end subroutine LU_mat

  subroutine LU_bsub_mat(self,b,trans)
    class(mat), intent(in) :: self
    type(mat), intent(inout) :: b
    logical, intent(in), optional :: trans
    character(len=1) :: trans_char
    integer(i4b) :: m,nrhs,info
    trans_char = 'N'    
    if(present(trans)) then
       if(trans) trans_char = 'T'
    end if
    m = self%m
    nrhs = b%n
    call check(self%square,'LU_bsub_mat','matrix not square')
    call check(self%real .eqv. b%real,'LU_bsub_mat','matrix and rhs must be of the same type')
    call check(self%factor == 'LU','LU_bsub_mat','matrix not LU factorised')
    call check(m == b%m,'LU_bsub_mat','rhs has the wrong dimensions')
    if(self%real) then
       call dgetrs(trans_char,m,nrhs,self%elem_r,m,self%ipiv,b%elem_r,m,info)
    else
       call zgetrs(trans_char,m,nrhs,self%elem_c,m,self%ipiv,b%elem_c,m,info)
    end if
    call check(info == 0,'LU_bsub_mat','problem with back substitution')
    return
  end subroutine LU_bsub_mat


  subroutine factor_mat(self,type)
    class(mat), intent(inout) :: self
    character(len=*), intent(in), optional :: type
    character(len=:), allocatable :: type_local
    if(present(type)) then
       type_local = type
    else
       type_local = 'LU'
    end if
    if(type_local == 'LU') then
       call self%LU()
    else
       call error('factor_mat','unknown scheme')
    end if
    return
  end subroutine factor_mat

  subroutine bsub_mat(self,b,trans)
    class(mat), intent(in) :: self
    type(mat), intent(inout) :: b
    logical, intent(in), optional :: trans
    call check(self%factor /= '','bsub_mat','matrix not factorised')
    if(self%factor == 'LU') then
       call self%LU_bsub(b,trans)
    else
       call error('bsub_mat','unknown factorisation scheme')
    end if
    return
  end subroutine bsub_mat
  


  !===============================================!
  !      procedures for band real matrices        !
  !===============================================!

  
  subroutine deallocate_bmat(self)
    class(bmat), intent(inout) :: self
    if(.not.self%allocated) return
    self%kl = 0
    self%ku = 0
    self%ld = 0
    self%band_set = .false.
    call self%mat%deallocate()
    return
  end subroutine deallocate_bmat

  subroutine band_bmat(self,kl,ku)
    class(bmat), intent(inout) :: self
    integer(i4b), intent(in) :: kl,ku
    call check(kl >= 0,'band_bmat','invalid lower bandwidth')
    call check(ku >= 0,'band_bmat','invalid upper bandwidth')
    self%kl = kl
    self%ku = ku
    self%band_set = .true.
    return
  end subroutine band_bmat
  

  subroutine allocate_bmat(self,m,n,real)
    class(bmat), intent(inout) :: self
    integer(i4b), intent(in) :: m,n
    logical, intent(in), optional :: real
    call check(self%band_set,'allocate_bmat','bandwidths not set')
    if(present(real)) self%real = real
    self%m  = m
    self%n  = n
    if(self%kl > m-1) self%kl = m-1
    if(self%ku > n-1) self%kl = n-1
    self%ld = 2*self%kl + self%ku + 1
    if(self%real) then       
       allocate(self%elem_r(self%ld,n))
       self%elem_r = 0.0_dp
    else
       allocate(self%elem_c(self%ld,n))
       self%elem_c = 0.0_dp
    end if
    if(m == n) self%square = .true.
    self%allocated = .true.
    return
  end subroutine allocate_bmat


  subroutine set_r_bmat(self,i,j,val)
    class(bmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'set_bmat','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'set_bmat','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'set_bmat','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    if(self%real) then
       self%elem_r(k,j) =  val
    else
       self%elem_c(k,j) =  val
    end if
    return
  end subroutine set_r_bmat


  subroutine set_c_bmat(self,i,j,val)
    class(bmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'set_bmat','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'set_bmat','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'set_bmat','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    if(self%real) then
       self%elem_r(k,j) =  val
    else
       self%elem_c(k,j) =  val
    end if
    return
  end subroutine set_c_bmat


  subroutine get_r_bmat(self,i,j,val)
    class(bmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(out) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'get_bmat','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'get_bmat','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'get_bmat','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    if(self%real) then
       val = self%elem_r(k,j)
    else
      val = self%elem_c(k,j) 
    end if
    return
  end subroutine get_r_bmat


  subroutine get_c_bmat(self,i,j,val)
    class(bmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(out) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'get_bmat','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'get_bmat','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'get_bmat','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    if(self%real) then
       val = self%elem_r(k,j) 
    else
       val = self%elem_c(k,j) 
    end if
    return
  end subroutine get_c_bmat

  
  subroutine inc_r_bmat(self,i,j,val)
    class(bmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'inc_bmat','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'inc_bmat','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'inc_bmat','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    if(self%real) then
       self%elem_r(k,j) = self%elem_r(k,j) + val
    else
       self%elem_c(k,j) = self%elem_c(k,j) + val
    end if
    return
  end subroutine inc_r_bmat


  subroutine inc_c_bmat(self,i,j,val)
    class(bmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'inc_bmat','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'inc_bmat','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'inc_bmat','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    if(self%real) then
       self%elem_r(k,j) = self%elem_r(k,j) + val
    else
       self%elem_c(k,j) = self%elem_c(k,j) + val
    end if
    return
  end subroutine inc_c_bmat

    
  subroutine LU_bmat(self)
    class(bmat), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'LU_bmat','matrix not allocated')
    if(self%factor == 'LU') return
    call check(self%factor == '','LU_mat','matrix already factored using different scheme')        
    allocate(self%ipiv(min(self%m,self%n)))
    if(self%real) then
       call dgbtrf(self%m,self%n,self%kl,self%ku,self%elem_r,self%ld,self%ipiv,info)
    else
       call zgbtrf(self%m,self%n,self%kl,self%ku,self%elem_c,self%ld,self%ipiv,info)
    end if
    call check(info == 0,'LU_bmat','problem with factorisation')
    self%factor = 'LU'
    return
  end subroutine LU_bmat
  

  subroutine LU_bsub_bmat(self,b,trans)
    class(bmat), intent(in) :: self
    type(mat), intent(inout) :: b
    logical, intent(in), optional :: trans
    character(len=1) :: trans_char
    integer(i4b) :: m,nrhs,info
    trans_char = 'N'    
    if(present(trans)) then
       if(trans) trans_char = 'T'
    end if
    m = self%m
    nrhs = b%n
    call check(self%square,'LU_bsub_bmat','matrix not square')
    call check(self%factor == 'LU','LU_bsub_mat','matrix not LU factorised')
    call check(self%real .eqv. b%real,'LU_bsub_mat','matrix and rhs must be of the same type')
    call check(m == b%m,'LU_bsub_bmat','rhs has the wrong dimensions')
    if(self%real) then
       call dgbtrs(trans_char,m,self%kl,self%ku,nrhs,self%elem_r,self%ld,self%ipiv,b%elem_r,m,info)
    else
       call zgbtrs(trans_char,m,self%kl,self%ku,nrhs,self%elem_c,self%ld,self%ipiv,b%elem_c,m,info)
    end if
    call check(info == 0,'LU_bsub_bmat','problem with back substitution')
    return
  end subroutine LU_bsub_bmat


  integer(i4b) function row_ind_bmat(self,i,j) result(k)
    class(bmat), intent(in) :: self
    integer(i4b), intent(in) :: i,j
    if(self%check) then
       call check(self%allocated,'inc_bmat','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'inc_bmat','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'inc_bmat','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j    
    return
  end function row_ind_bmat



  !============================================================!
  !    procedures for positive symmetric band real matrices    !
  !============================================================!
    
  subroutine deallocate_hbmat(self)
    class(hbmat), intent(inout) :: self
    if(.not.self%allocated) return
    self%kd = 0
    self%ld = 0
    call self%mat%deallocate()
    return
  end subroutine deallocate_hbmat


  subroutine band_hbmat(self,kd)
    class(hbmat), intent(inout) :: self
    integer(i4b), intent(in) :: kd
    call check(kd >= 0,'band_hbmat','invalid bandwidth')
    self%kd = kd
    self%band_set = .true.
    return
  end subroutine band_hbmat
  

  subroutine allocate_hbmat(self,m,n,real)
    class(hbmat), intent(inout) :: self
    integer(i4b), intent(in) :: m,n
    logical, intent(in), optional :: real
    call check(m == n,'allocate_hbmat','matices must be square')
    call check(self%band_set,'allocate_hbmat','bandwidths not set')
    if(present(real)) self%real = real
    self%m  = m
    self%n  = n
    if(self%kd > m-1) self%kd = m-1
    self%ld = self%kd+1
    if(self%real) then
       allocate(self%elem_r(self%ld,m))
       self%elem_r = 0.0_dp
    else
       allocate(self%elem_c(self%ld,m))
       self%elem_c = 0.0_dp
    end if
    self%square = .true.
    self%allocated = .true.
    return
  end subroutine allocate_hbmat

    
  subroutine set_r_hbmat(self,i,j,val)
    class(hbmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    integer(i4b) :: k,il,jl
    if(il <= jl) then
       il = i
       jl = j
    else
       il = j
       jl = i
    end if
    if(self%check) then
       call check(self%allocated,'set_hbmat','matrix not allocated')
       call check(il >= max(1,jl-self%kd),'set_hbmat','row index out of range')
       call check(jl >= 1 .and.  jl <= self%m,'set_hbmat','column index out of range')

    end if
    k = self%kd + 1 + i - j
    if(self%real) then
       self%elem_r(k,j) =  val
    else
       self%elem_c(k,j) =  val
    end if
    return
  end subroutine set_r_hbmat


  subroutine set_c_hbmat(self,i,j,val)
    class(hbmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val
    integer(i4b) :: k,il,jl
    if(il <= jl) then
       il = i
       jl = j
    else
       il = j
       jl = i
    end if
    if(self%check) then
       call check(self%allocated,'set_hbmat','matrix not allocated')
       call check(il >= max(1,jl-self%kd),'set_hbmat','row index out of range')
       call check(jl >= 1 .and.  jl <= self%m,'set_hbmat','column index out of range')

    end if    
    k = self%kd + 1 + i - j
    if(self%real) then
       self%elem_r(k,j) =  val
    else
       self%elem_c(k,j) =  val
    end if
    return
  end subroutine set_c_hbmat


  subroutine get_r_hbmat(self,i,j,val)
    class(hbmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(out) :: val
    integer(i4b) :: k,il,jl
    if(il <= jl) then
       il = i
       jl = j
    else
       il = j
       jl = i
    end if
    if(self%check) then
       call check(self%allocated,'get_hbmat','matrix not allocated')
       call check(il >= max(1,jl-self%kd),'get_hbmat','row index out of range')
       call check(jl >= 1 .and.  jl <= self%m,'get_hbmat','column index out of range')

    end if
    k = self%kd + 1 + i - j
    if(self%real) then
       val = self%elem_r(k,j) 
    else
       val = self%elem_c(k,j)
    end if
    return
  end subroutine get_r_hbmat


  subroutine get_c_hbmat(self,i,j,val)
    class(hbmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(out) :: val
    integer(i4b) :: k,il,jl
    if(il <= jl) then
       il = i
       jl = j
    else
       il = j
       jl = i
    end if
    if(self%check) then
       call check(self%allocated,'get_hbmat','matrix not allocated')
       call check(il >= max(1,jl-self%kd),'get_hbmat','row index out of range')
       call check(jl >= 1 .and.  jl <= self%m,'get_hbmat','column index out of range')

    end if    
    k = self%kd + 1 + i - j
    if(self%real) then
       val = self%elem_r(k,j) 
    else
       val = self%elem_c(k,j) 
    end if
    return
  end subroutine get_c_hbmat
  
  subroutine inc_r_hbmat(self,i,j,val)
    class(hbmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    integer(i4b) :: k,il,jl
    if(il <= jl) then
       il = i
       jl = j
    else
       il = j
       jl = i
    end if
    if(self%check) then
       call check(self%allocated,'set_hbmat','matrix not allocated')
       call check(il >= max(1,jl-self%kd),'set_hbmat','row index out of range')
       call check(jl >= 1 .and.  jl <= self%m,'set_hbmat','column index out of range')

    end if    
    k = self%kd + 1 + i - j
    if(self%real) then
       self%elem_r(k,j) =  self%elem_r(k,j) + val
    else
       self%elem_c(k,j) =  self%elem_c(k,j) + val
    end if
    return
  end subroutine inc_r_hbmat


  subroutine inc_c_hbmat(self,i,j,val)
    class(hbmat), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: val
    integer(i4b) :: k,il,jl
    if(il <= jl) then
       il = i
       jl = j
    else
       il = j
       jl = i
    end if
    if(self%check) then
       call check(self%allocated,'set_hbmat','matrix not allocated')
       call check(il >= max(1,jl-self%kd),'set_hbmat','row index out of range')
       call check(jl >= 1 .and.  jl <= self%m,'set_hbmat','column index out of range')

    end if    
    k = self%kd + 1 + i - j
    if(self%real) then
       self%elem_r(k,j) =  self%elem_r(k,j) + val
    else
       self%elem_c(k,j) =  self%elem_c(k,j) + val
    end if
    return
  end subroutine inc_c_hbmat

  
  subroutine LU_hbmat(self)
    class(hbmat), intent(inout) :: self
    call error('LU_hbmat','factorisation not implemented')
    return
  end subroutine LU_hbmat

  subroutine LU_bsub_hbmat(self,b,trans)
    class(hbmat), intent(in) :: self
    type(mat), intent(inout) :: b
    logical, intent(in), optional :: trans
    call error('LU_bsub_hbmat','factorisation not implemented')
    return
  end subroutine LU_bsub_hbmat


  integer(i4b) function row_ind_hbmat(self,i,j) result(k)
    class(hbmat), intent(in) :: self
    integer(i4b), intent(in) :: i,j
    integer(i4b) :: il,jl
    if(il <= jl) then
       il = i
       jl = j
    else
       il = j
       jl = i
    end if
    if(self%check) then
       call check(self%allocated,'set_hbmat','matrix not allocated')
       call check(il >= max(1,jl-self%kd),'set_hbmat','row index out of range')
       call check(jl >= 1 .and.  jl <= self%m,'set_hbmat','column index out of range')

    end if    
    k = self%kd + 1 + i - j
    return
  end function row_ind_hbmat

  
  subroutine cholesky_hbmat(self)
    class(hbmat), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'cholesky_bmat','matrix not allocated')
    if(self%factor == 'cholesky') return
    call check(self%factor == '','cholesky_hbmat','matrix already factored using different scheme')
    if(self%real) then
       call dpbtrf('U',self%n,self%kd,self%elem_r,self%ld,info)
    else
       call zpbtrf('U',self%n,self%kd,self%elem_c,self%ld,info)
    end if
    call check(info == 0,'cholesky_bmat','problem with factorisation')
    self%factor = 'cholesky'
    return
  end subroutine cholesky_hbmat


  subroutine cholesky_bsub_hbmat(self,b,trans)
    class(hbmat), intent(in) :: self
    type(mat), intent(inout) :: b
    logical, intent(in), optional :: trans
    integer(i4b) :: info
    call check(self%allocated,'cholesky_bsub_bmat','matrix not allocated')
    call check(self%factor == 'cholesky','cholesky_bsub_bmat','matrix not factorised')
    call check(self%real .eqv. b%real,'cholesky_bsub_mat','matrix and rhs must be of the same type')
    call check(self%m == b%m,'cholesky_bsub_bmat','rhs has the wrong dimensions')
    if(self%real) then
       call dpbtrs('U',self%n,self%kd,b%n,self%elem_r,self%ld,b%elem_r,b%m,info)
    else
       call zpbtrs('U',self%n,self%kd,b%n,self%elem_c,self%ld,b%elem_c,b%m,info)
    end if
    call check(info == 0,'cholesky_bmat','problem with back substitution')
    return
  end subroutine Cholesky_bsub_hbmat

 
  subroutine factor_hbmat(self,type)
    class(hbmat), intent(inout) :: self
    character(len=*), intent(in), optional :: type
    character(len=:), allocatable :: type_local
    if(present(type)) then
       type_local = type
    else
       type_local = 'cholesky'
    end if
    if(type_local == 'cholesky') then
       call self%Cholesky()
    else if(type_local == 'LU') then
       call self%LU()
    else
       call error('factor_hbmat','unknown scheme')
    end if
    return
  end subroutine factor_hbmat

  subroutine bsub_hbmat(self,b,trans)
    class(hbmat), intent(in) :: self
    type(mat), intent(inout) :: b
    logical, intent(in), optional :: trans
    call check(self%factor /= '','bsub_hbmat','matrix not factorised')
    if(self%factor == 'LU') then
       call self%LU_bsub(b,trans)
    else if(self%factor == 'cholesky') then
       call self%Cholesky_bsub(b,trans)
    else
       call error('bsub_hbmat','unknown factorisation scheme')
    end if
    return
  end subroutine bsub_hbmat


  
  
end module module_LAPACK
