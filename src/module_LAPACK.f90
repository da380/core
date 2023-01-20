 module module_LAPACK

  use module_constants
  use module_error
  implicit none

  type rm
     logical :: allocated = .false.
     logical :: square = .false.
     logical :: check = .false.
     character(len=:), allocatable, private :: factor  
     integer(i4b) :: m = 0
     integer(i4b) :: n = 0
     real(dp), dimension(:,:), allocatable :: elem
     integer(i4b), dimension(:), allocatable, private :: ipiv
   contains
     procedure :: deallocate => deallocate_rm
     procedure :: allocate => allocate_rm
     procedure :: set        => set_rm
     procedure :: inc        => inc_rm
     procedure :: set_sym    => set_sym_rm
     procedure :: inc_sym    => inc_sym_rm
     procedure :: LU         => LU_rm
     procedure :: LU_bsub => LU_bsub_rm
     procedure :: fac        => factor_rm
     procedure :: bsub       => bsub_rm
  end type rm

  type, extends(rm) :: brm
     logical :: band_set = .false.
     integer(i4b) :: kl
     integer(i4b) :: ku
     integer(i4b) :: ld
   contains
     procedure :: deallocate => deallocate_brm
     procedure :: band => band_brm
     procedure :: allocate => allocate_brm
     procedure :: set        => set_brm
     procedure :: inc        => inc_brm
     procedure :: set_sym    => set_sym_brm
     procedure :: inc_sym    => inc_sym_brm
     procedure :: LU         => LU_brm
     procedure :: LU_bsub => LU_bsub_brm
     procedure :: row_ind    => row_ind_brm
  end type brm

  type, extends(rm) :: psbrm
     logical :: band_set = .false.
     integer(i4b) :: kd
     integer(i4b) :: ld
   contains
     procedure :: deallocate     => deallocate_psbrm
     procedure :: band => band_psbrm
     procedure :: allocate => allocate_psbrm
     procedure :: set            => set_psbrm
     procedure :: inc            => inc_psbrm
     procedure :: set_sym        => set_psbrm
     procedure :: inc_sym        => inc_psbrm
     procedure :: LU             => LU_psbrm
     procedure :: LU_bsub        => LU_bsub_psbrm
     procedure :: row_ind        => row_ind_psbrm
     procedure :: Cholesky       => Cholesky_psbrm
     procedure :: Cholesky_bsub  => Cholesky_bsub_psbrm
     procedure :: fac            => factor_psbrm
     procedure :: bsub           => bsub_psbrm
  end type psbrm


  
contains


  !===============================================!
  !          procedures for real matrices         !
  !===============================================!
    
  subroutine deallocate_rm(self)
    class(rm), intent(inout) :: self
    if(.not.self%allocated) return
    self%m = 0
    self%n = 0
    if(allocated(self%ipiv)) deallocate(self%ipiv)
    if(allocated(self%elem)) deallocate(self%elem)
    if(allocated(self%factor)) deallocate(self%factor)
    self%allocated = .false.    
    return
  end subroutine deallocate_rm


  subroutine allocate_rm(self,m,n)
    class(rm), intent(inout) :: self
    integer(i4b), intent(in) :: m,n
    self%m = m
    self%n = n
    allocate(self%elem(m,n))
    self%elem = 0.0_dp
    if(m == n) self%square = .true.
    self%allocated = .true.
    return
  end subroutine allocate_rm

  
  subroutine set_rm(self,i,j,val)
    class(rm), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'set_rm','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'set_rm','row index out of range')
       call check(j >= 1 .and. j <= self%n,'set_rm','column index out of range')
    end if
    self%elem(i,j) =  val
    return
  end subroutine set_rm
  
  subroutine inc_rm(self,i,j,val)
    class(rm), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'inc_rm','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'inc_rm','row index out of range')
       call check(j >= 1 .and. j <= self%n,'inc_rm','column index out of range')
    end if
    self%elem(i,j) = self%elem(i,j) + val
    return
  end subroutine inc_rm


  subroutine set_sym_rm(self,i,j,val)
    class(rm), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'set_rm','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'set_sym_rm','row index out of range')
       call check(j >= 1 .and. j <= self%n,'set_sym_rm','column index out of range')
    end if    
    self%elem(i,j) = val
    if(i /= j) then
       self%elem(j,i) =  val
    end if
    return
  end subroutine set_sym_rm

  
  
  subroutine inc_sym_rm(self,i,j,val)
    class(rm), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    if(self%check) then
       call check(self%allocated,'inc_rm','matrix not allocated')
       call check(i >= 1 .and. i <= self%m,'inc_rm','row index out of range')
       call check(j >= 1 .and. j <= self%n,'inc_rm','column index out of range')
    end if    
    self%elem(i,j) = self%elem(i,j) + val
    if(i /= j) then
       self%elem(j,i) = self%elem(j,i) + val
    end if
    return
  end subroutine inc_sym_rm
  
    
  subroutine LU_rm(self)
    class(rm), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'LU_rm','matrix not allocated')
    if(self%factor == 'LU') return
    call check(self%factor == '','LU_rm','matrix already factored using different scheme')    
    allocate(self%ipiv(min(self%m,self%n)))
    call dgetrf(self%m,self%n,self%elem,self%m,self%ipiv,info)
    call check(info == 0,'LU_rm','problem with factorisation')
    self%factor = 'LU'
    return
  end subroutine LU_rm

  subroutine LU_bsub_rm(self,b,trans)
    class(rm), intent(in) :: self
    type(rm), intent(inout) :: b
    logical, intent(in), optional :: trans
    character(len=1) :: trans_char
    integer(i4b) :: m,nrhs,info
    trans_char = 'N'    
    if(present(trans)) then
       if(trans) trans_char = 'T'
    end if
    m = self%m
    nrhs = b%n
    call check(self%square,'LU_bsub_rm','matrix not square')
    call check(self%factor == 'LU','LU_bsub_rm','matrix not LU factorised')
    call check(m == b%m,'LU_bsub_rm','rhs has the wrong dimensions')
    call dgetrs	(trans_char,m,nrhs,self%elem,m,self%ipiv,b%elem,m,info)
    call check(info == 0,'LU_bsub_rm','problem with back substitution')
    return
  end subroutine LU_bsub_rm


  subroutine factor_rm(self,type)
    class(rm), intent(inout) :: self
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
       call error('factor_rm','unknown scheme')
    end if
    return
  end subroutine factor_rm

  subroutine bsub_rm(self,b,trans)
    class(rm), intent(in) :: self
    type(rm), intent(inout) :: b
    logical, intent(in), optional :: trans
    call check(self%factor /= '','bsub_rm','matrix not factorised')
    if(self%factor == 'LU') then
       call self%LU_bsub(b,trans)
    else
       call error('bsub_rm','unknown factorisation scheme')
    end if
    return
  end subroutine bsub_rm
  


  !===============================================!
  !      procedures for band real matrices        !
  !===============================================!
    
  subroutine deallocate_brm(self)
    class(brm), intent(inout) :: self
    if(.not.self%allocated) return
    self%kl = 0
    self%ku = 0
    self%ld = 0
    self%band_set = .false.
    call self%rm%deallocate()
    return
  end subroutine deallocate_brm

  subroutine band_brm(self,kl,ku)
    class(brm), intent(inout) :: self
    integer(i4b), intent(in) :: kl,ku
    call check(kl >= 0,'band_brm','invalid lower bandwidth')
    call check(ku >= 0,'band_brm','invalid upper bandwidth')
    self%kl = kl
    self%ku = ku
    self%band_set = .true.
    return
  end subroutine band_brm
  

  subroutine allocate_brm(self,m,n)
    class(brm), intent(inout) :: self
    integer(i4b), intent(in) :: m,n
    call check(self%band_set,'allocate_brm','bandwidths not set')
    self%m  = m
    self%n  = n
    if(self%kl > m-1) self%kl = m-1
    if(self%ku > n-1) self%kl = n-1
    self%ld = 2*self%kl + self%ku + 1
    allocate(self%elem(self%ld,n))
    self%elem = 0.0_dp
    if(m == n) self%square = .true.
    self%allocated = .true.
    return
  end subroutine allocate_brm


  subroutine set_brm(self,i,j,val)
    class(brm), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'set_brm','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'set_brm','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'set_brm','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    self%elem(k,j) =  val
    return
  end subroutine set_brm
  
  subroutine inc_brm(self,i,j,val)
    class(brm), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'inc_brm','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'inc_brm','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'inc_brm','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    self%elem(k,j) = self%elem(k,j) + val
    return
  end subroutine inc_brm


  subroutine set_sym_brm(self,i,j,val)
    class(brm), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'set_brm','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'set_sym_brm','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'set_sym_brm','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    self%elem(k,j) =  val
    if(i /= j) then
       k = self%kl + self%ku + 1 + j - i
       self%elem(k,i) =  val
    end if
    return
  end subroutine set_sym_brm

  subroutine inc_sym_brm(self,i,j,val)
    class(brm), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: val
    integer(i4b) :: k
    if(self%check) then
       call check(self%allocated,'inc_brm','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'inc_sym_brm','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'inc_sym_brm','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j
    self%elem(k,j) = self%elem(k,j) + val
    if(i /= j) then
       k = self%kl + self%ku + 1 + j - i
       self%elem(k,i) = self%elem(k,i) + val
    end if
    return
  end subroutine inc_sym_brm
    
  subroutine LU_brm(self)
    class(brm), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'LU_brm','matrix not allocated')
    if(self%factor == 'LU') return
    call check(self%factor == '','LU_rm','matrix already factored using different scheme')        
    allocate(self%ipiv(min(self%m,self%n)))
    call dgbtrf(self%m,self%n,self%kl,self%ku,self%elem,self%ld,self%ipiv,info)
    call check(info == 0,'LU_brm','problem with factorisation')
    self%factor = 'LU'
    return
  end subroutine LU_brm

  subroutine LU_bsub_brm(self,b,trans)
    class(brm), intent(in) :: self
    type(rm), intent(inout) :: b
    logical, intent(in), optional :: trans
    character(len=1) :: trans_char
    integer(i4b) :: m,nrhs,info
    trans_char = 'N'    
    if(present(trans)) then
       if(trans) trans_char = 'T'
    end if
    m = self%m
    nrhs = b%n
    call check(self%square,'LU_bsub_brm','matrix not square')
    call check(self%factor == 'LU','LU_bsub_rm','matrix not LU factorised')
    call check(m == b%m,'LU_bsub_brm','rhs has the wrong dimensions')
    call dgbtrs(trans_char,m,self%kl,self%ku,nrhs,self%elem,self%ld,self%ipiv,b%elem,m,info)	
    call check(info == 0,'LU_bsub_brm','problem with back substitution')
    return
  end subroutine LU_bsub_brm


  integer(i4b) function row_ind_brm(self,i,j) result(k)
    class(brm), intent(in) :: self
    integer(i4b), intent(in) :: i,j
    if(self%check) then
       call check(self%allocated,'inc_brm','matrix not allocated')
       call check(i >= 1 .and. i <= min(self%m,j+self%kl),'inc_brm','row index out of range')
       call check(j >= 1 .and. j <= min(self%n,i+self%kl),'inc_brm','column index out of range')       
    end if    
    k = self%kl + self%ku + 1 + i - j    
    return
  end function row_ind_brm



  !============================================================!
  !      procedures for positive symmetric band real matrices  !
  !============================================================!
    
  subroutine deallocate_psbrm(self)
    class(psbrm), intent(inout) :: self
    if(.not.self%allocated) return
    self%kd = 0
    self%ld = 0
    call self%rm%deallocate()
    return
  end subroutine deallocate_psbrm


  subroutine band_psbrm(self,kd)
    class(psbrm), intent(inout) :: self
    integer(i4b), intent(in) :: kd
    call check(kd >= 0,'band_psbrm','invalid bandwidth')
    self%kd = kd
    self%band_set = .true.
    return
  end subroutine band_psbrm
  

  subroutine allocate_psbrm(self,m,n)
    class(psbrm), intent(inout) :: self
    integer(i4b), intent(in) :: m,n
    call check(m == n,'allocate_psbrm','matices must be square')
    call check(self%band_set,'allocate_psbrm','bandwidths not set')
    self%m  = m
    self%n  = n
    if(self%kd > m-1) self%kd = m-1
    self%ld = self%kd+1
    allocate(self%elem(self%ld,m))
    self%elem = 0.0_dp
    self%square = .true.
    self%allocated = .true.
    return
  end subroutine allocate_psbrm

    
  subroutine set_psbrm(self,i,j,val)
    class(psbrm), intent(inout) :: self
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
       call check(self%allocated,'set_psbrm','matrix not allocated')
       call check(il >= max(1,jl-self%kd),'set_psbrm','row index out of range')
       call check(jl >= 1 .and.  jl <= self%m,'set_psbrm','column index out of range')

    end if    
    k = self%kd + 1 + i - j
    self%elem(k,j) =  val
    return
  end subroutine set_psbrm
  
  subroutine inc_psbrm(self,i,j,val)
    class(psbrm), intent(inout) :: self
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
       call check(self%allocated,'set_psbrm','matrix not allocated')
       call check(il >= max(1,jl-self%kd),'set_psbrm','row index out of range')
       call check(jl >= 1 .and.  jl <= self%m,'set_psbrm','column index out of range')

    end if    
    k = self%kd + 1 + i - j
    self%elem(k,j) =  self%elem(k,j) + val
    return
  end subroutine inc_psbrm
    
  subroutine LU_psbrm(self)
    class(psbrm), intent(inout) :: self
    call error('LU_psbrm','factorisation not implemented')
    return
  end subroutine LU_psbrm

  subroutine LU_bsub_psbrm(self,b,trans)
    class(psbrm), intent(in) :: self
    type(rm), intent(inout) :: b
    logical, intent(in), optional :: trans
    call error('LU_bsub_psbrm','factorisation not implemented')
    return
  end subroutine LU_bsub_psbrm


  integer(i4b) function row_ind_psbrm(self,i,j) result(k)
    class(psbrm), intent(in) :: self
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
       call check(self%allocated,'set_psbrm','matrix not allocated')
       call check(il >= max(1,jl-self%kd),'set_psbrm','row index out of range')
       call check(jl >= 1 .and.  jl <= self%m,'set_psbrm','column index out of range')

    end if    
    k = self%kd + 1 + i - j
    return
  end function row_ind_psbrm

  

  subroutine Cholesky_psbrm(self)
    class(psbrm), intent(inout) :: self
    integer(i4b) :: info
    call check(self%allocated,'Cholesky_brm','matrix not allocated')
    if(self%factor == 'Cholesky') return
    call check(self%factor == '','Cholesky_psbrm','matrix already factored using different scheme')        
    call dpbtrf('U',self%n,self%kd,self%elem,self%ld,info)
    call check(info == 0,'Cholesky_brm','problem with factorisation')
    self%factor = 'Cholesky'
    return
  end subroutine Cholesky_psbrm


  subroutine Cholesky_bsub_psbrm(self,b,trans)
    class(psbrm), intent(in) :: self
    type(rm), intent(inout) :: b
    logical, intent(in), optional :: trans
    integer(i4b) :: info
    call check(self%allocated,'Cholesky_bsub_brm','matrix not allocated')
    call check(self%factor == 'Cholesky','Cholesky_bsub_brm','matrix not factorised')
    call check(self%m == b%m,'Cholesky_bsub_brm','rhs has the wrong dimensions')
    call dpbtrs	('U',self%n,self%kd,b%n,self%elem,self%ld,b%elem,b%m,info)
    call check(info == 0,'Cholesky_brm','problem with back substitution')
    return
  end subroutine Cholesky_bsub_psbrm

 
  subroutine factor_psbrm(self,type)
    class(psbrm), intent(inout) :: self
    character(len=*), intent(in), optional :: type
    character(len=:), allocatable :: type_local
    if(present(type)) then
       type_local = type
    else
       type_local = 'Cholesky'
    end if
    if(type_local == 'Cholesky') then
       call self%Cholesky()
    else if(type_local == 'LU') then
       call self%LU()
    else
       call error('factor_psbrm','unknown scheme')
    end if
    return
  end subroutine factor_psbrm

  subroutine bsub_psbrm(self,b,trans)
    class(psbrm), intent(in) :: self
    type(rm), intent(inout) :: b
    logical, intent(in), optional :: trans
    call check(self%factor /= '','bsub_psbrm','matrix not factorised')
    if(self%factor == 'LU') then
       call self%LU_bsub(b,trans)
    else if(self%factor == 'Cholesky') then
       call self%Cholesky_bsub(b,trans)
    else
       call error('bsub_psbrm','unknown factorisation scheme')
    end if
    return
  end subroutine bsub_psbrm
  
  
end module module_LAPACK
