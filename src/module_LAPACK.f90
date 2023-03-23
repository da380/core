module module_LAPACK

  use module_constants
  use module_error
  implicit none


  
  
  type MatLAP
     !---------------------------------------------------------------------!
     !---------------------------------------------------------------------!
     ! dr  = dense, real
     ! dc  = dense, complex
     ! dsr = dense, symmetric, real    
     ! dsc = dense, symmetric, complex 
     ! dsp = dense, symmetric, real, positive-definite
     ! dhp = dense, hermitian, complex, positive-definite
     ! br  = banded, real
     ! bc  = banded, complex
     ! bsp = banded, symmetric, real, positive-definite
     ! bhp = banded, hermitian, complex, positive-definite
     !---------------------------------------------------------------------!
     !---------------------------------------------------------------------!
     private
     logical :: allocated = .false.
     logical :: factored = .false.
     character(len=:), allocatable :: type
     integer(i4b) :: m
     integer(i4b) :: n
     integer(i4b) :: kl 
     integer(i4b) :: ku
     integer(i4b) :: kd
     integer(i4b) :: ld
     real(dp), dimension(:,:), allocatable :: data_r
     complex(dpc), dimension(:,:), allocatable :: data_c
     integer(i4b), dimension(:), allocatable, private :: ipiv
   contains
     procedure :: deallocate => deallocate_MatLAP
     procedure :: allocate   => allocate_MatLAP
     procedure :: row_dim => row_dim_MatLAP
     procedure :: col_dim => col_dim_MatLAP
     procedure :: square => square_MatLAP
     procedure :: mirror_upper   => mirror_upper_MatLAP
     procedure :: set_value_MatLAP_r
     procedure :: set_value_MatLAP_c
     procedure :: set_value_MatLAP_single_r
     procedure :: set_value_MatLAP_single_c
     generic   :: set => set_value_MatLAP_r,        &
                         set_value_MatLAP_single_r, &
                         set_value_MatLAP_c,        &
                         set_value_MatLAP_single_c
     procedure :: factor     => factor_MatLAP
     procedure, private :: backsub_MatLAP_r_single
     procedure, private :: backsub_MatLAP_r_many
     procedure, private :: backsub_MatLAP_c_single
     procedure, private :: backsub_MatLAP_c_many
     generic   :: backsub => backsub_MatLAP_r_single, &
                             backsub_MatLAP_r_many,   &
                             backsub_MatLAP_c_single, &
                             backsub_MatLAP_c_many
  end type MatLAP

  
  
contains


  
  !=====================================================================!
  !                      routines for MatLAP type                       !
  !=====================================================================!



  subroutine deallocate_MatLAP(self)
    class(MatLAP), intent(inout) :: self
    if(.not.self%allocated)  return
    if(allocated(self%type)) deallocate(self%type)
    if(allocated(self%data_r)) deallocate(self%data_r)
    if(allocated(self%data_c)) deallocate(self%data_c)
    if(allocated(self%ipiv)) deallocate(self%ipiv)
    self%allocated = .false.
    return
  end subroutine deallocate_MatLAP

  subroutine allocate_MatLAP(self,type,m,n,kl,ku,kd)
    class(MatLAP), intent(inout) :: self
    character(len=*), intent(in) :: type
    integer(i4b), intent(in) :: m
    integer(i4b), intent(in) :: n
    integer(i4b), intent(in), optional :: kl
    integer(i4b), intent(in), optional :: ku
    integer(i4b), intent(in), optional :: kd
    
    logical :: test
    
    call self%deallocate()    
    call check(m > 0,'allocate_MatLAP','non-positive row dimension')
    call check(n > 0,'allocate_MatLAP','non-positive column dimension')
    self%type = type
    self%m = m
    self%n = n
    
    if(type == 'dr') then
       
       allocate(self%data_r(m,n))
       self%data_r(m,n) = 0.0_dp
       
    else if(type == 'dsr' .or. type == 'dsp') then
       
       call check(m == n,'allocate_MatLAP','symmetric/hermitian matrices must be square')
       allocate(self%data_r(m,n))
       self%data_r(m,n) = 0.0_dp
       
    else if(type == 'dc') then
       
       allocate(self%data_c(m,n))
       self%data_c(m,n) = 0.0_dp
       
    else if(type == 'dsc' .or. type == 'dh' .or. type == 'dhp') then
       
       call check(m == n,'allocate_MatLAP','symmetric/hermitian matrices must be square')
       allocate(self%data_c(m,n))
       self%data_c(m,n) = 0.0_dp
       
    else if(type == 'br' .or. type == 'bc') then
       
       test = present(kl) .and. present(ku)
       call check(test,'allocate_MatLAP','missing information for banded matrices')
       call check(kl >= 0,'allocate_MatLAP','negative lower bandwidth')
       call check(ku >= 0,'allocate_MatLAP','negative upper bandwidth')       
       self%kl = kl
       self%ku = ku
       self%ld = 2*kl+ku+1
       if(type == 'br') then
          allocate(self%data_r(self%ld,n))
          self%data_r = 0.0_dp       
       else
          allocate(self%data_c(self%ld,n))
          self%data_c = 0.0_dp       
       end if
       
    else if(type == 'bsp' .or. type == 'bhp') then
       
       call check(m == n,'allocate_MatLAP','symmetric/hermitian matrices must be square')
       test = present(kd)
       call check(test,'allocate_MatLAP','missing information for hermitian banded matrices')
       call check(kd >= 0,'allocate_MatLAP','negative bandwidth')
       self%kd = kd
       self%ld = kd+1
       if(type == 'bsp') then
          allocate(self%data_r(self%ld,n))
          self%data_r = 0.0_dp       
       else
          allocate(self%data_c(self%ld,n))
          self%data_c = 0.0_dp       
       end if
       
    else
       
       call error('allocate_MatLAP','undefined type')
       
    end if

    self%allocated = .true.
    
    return
  end subroutine allocate_MatLAP


  integer(i4b) function row_dim_MatLAP(self) result(m)
    class(MatLAP), intent(in) :: self
    m = self%m
    return
  end function row_dim_MatLAP


  integer(i4b) function col_dim_MatLAP(self) result(n)
    class(MatLAP), intent(in) :: self
    n = self%n
    return
  end function col_dim_MatLAP

  
  logical function square_MatLAP(self) result(bool)
    class(MatLAP), intent(in) :: self
    bool = self%m == self%n
    return
  end function square_MatLAP


  subroutine mirror_upper_MatLAP(self)
    class(MatLAP), intent(inout) :: self

    integer(i4b) :: i,j,k1,k2

    call check(self%allocated,'mirror_upper_MatLAP','matrix not allocated')
    
    associate(type => self%type,    &
              m    => self%m,       &
              n    => self%n,       &
              kl   => self%kl,      &
              ku   => self%ku,      &
              kd   => self%kd,      &
              ld   => self%ld,      &              
              ar   => self%data_r,  &
              ac   => self%data_c)

      if(type == 'dr') then

         do j = 1,n
            do i = j+1,m
               ar(i,j) = ar(j,i)
            end do
         end do

      else if(type == 'dc') then

         do j = 1,n
            do i = j+1,m
               ac(i,j) = ac(j,i)
            end do
         end do

      else if(type == 'br') then

         do j = 1,n
            do i = j+1,min(m,j+kl)
               k1 = kl + ku + 1 + i - j
               k2 = kl + ku + 1 + j - i
               ar(k1,j) = ar(k2,i)
            end do
         end do

      else if(type == 'bc') then

         do j = 1,n
            do i = j+1,min(m,j+kl)
               k1 = kl + ku + 1 + i - j
               k2 = kl + ku + 1 + j - i
               ac(k1,j) = ac(k2,i)
            end do
         end do         
         
      end if
      
    end associate

    
    return
  end subroutine mirror_upper_MatLAP
  
  subroutine set_value_MatLAP_r(self,ia,ma,ja,na,a,inc)
    class(MatLAP), intent(inout) :: self
    integer(i4b), intent(in) :: ia,ma,ja,na
    real(dp), dimension(ma,na), intent(in) :: a
    logical, intent(in), optional  :: inc

    logical :: incl
    integer(i4b) :: i,j,k


    if(present(inc)) then
       incl = inc
    else
       incl = .false.
    end if
    
    call check(self%allocated,'set_value_MatLAP_r','matrix not allocated')
    
    associate(type => self%type,    &
              m    => self%m,       &
              n    => self%n,       &
              kl   => self%kl,      &
              ku   => self%ku,      &
              kd   => self%kd,      &
              ld   => self%ld,      &              
              ar   => self%data_r,  &
              ac   => self%data_c)

      call check(ia >= 1 .and. ia+ma-1 <= m,'set_value_MatLAP_r','row indices out of range')
      call check(ja >= 1 .and. ja+na-1 <= n,'set_value_MatLAP_r','column indices out of range')

      if(type == 'dr') then
         
         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               if(incl) then
                  ar(i,j) = ar(i,j) + a(i-ia+1,j-ja+1)
               else
                  ar(i,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do
         
      else if( type == 'dsr' .or. type == 'dsp') then
         
         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               if(i > j) cycle
               if(incl) then
                  ar(i,j) = ar(i,j) + a(i-ia+1,j-ja+1)
               else
                  ar(i,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do
         
      else if(type == 'dc') then
         
         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               if(incl) then
                  ac(i,j) = ac(i,j) + a(i-ia+1,j-ja+1)
               else
                  ac(i,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do

      else if(type == 'dsc' .or. type == 'dhp') then

         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               if(i > j) cycle
               if(incl) then
                  ac(i,j) = ac(i,j) + a(i-ia+1,j-ja+1)
               else
                  ac(i,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do

      else if(type == 'br') then

         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               k = kl + ku + 1 + i - j
               if(k < 1 .or. k > ld) cycle
               if(incl) then
                  ar(k,j) = ar(k,j) + a(i-ia+1,j-ja+1)
               else
                  ar(k,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do

      else if(type == 'bc') then

         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               k = kl + ku + 1 + i - j
               if(k < 1 .or. k > ld) cycle
               if(incl) then
                  ac(k,j) = ac(k,j) + a(i-ia+1,j-ja+1)
               else
                  ac(k,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do         

      else if(type == 'bsp') then
         
         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               k = kd + 1 + i - j
               if(k < 1 .or. k > ld) cycle
               if(incl) then
                  ar(k,j) = ar(k,j) + a(i-ia+1,j-ja+1)
               else
                  ar(k,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do

      else if(type == 'bhp') then
         
         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               k = kd + 1 + i - j
               if(k < 1 .or. k > ld) cycle
               if(incl) then
                  ac(k,j) = ac(k,j) + a(i-ia+1,j-ja+1)
               else
                  ac(k,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do
         
      else

         call error('set_value_MatLAP_r','unknown matrix type')
         
      end if
      
    end associate
         
    return
  end subroutine set_value_MatLAP_r


  subroutine set_value_MatLAP_single_r(self,i,j,a,inc)
    class(MatLAP), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    real(dp), intent(in) :: a
    logical, intent(in), optional :: inc
    real(dp), dimension(1,1) :: al
    al = a
    call self%set(i,1,j,1,al,inc)
    return
  end subroutine set_value_MatLAP_single_r

  
  subroutine set_value_MatLAP_c(self,ia,ma,ja,na,a,inc)
    class(MatLAP), intent(inout) :: self
    integer(i4b), intent(in) :: ia,ma,ja,na
    complex(dpc), dimension(ma,na), intent(in) :: a
    logical, intent(in), optional  :: inc

    logical :: incl
    integer(i4b) :: i,j,k


    if(present(inc)) then
       incl = inc
    else
       incl = .false.
    end if
    
    call check(self%allocated,'set_value_MatLAP_c','matrix not allocated')
    
    associate(type => self%type,    &
              m    => self%m,       &
              n    => self%n,       &
              kl   => self%kl,      &
              ku   => self%ku,      &
              kd   => self%kd,      &
              ld   => self%ld,      &              
              ac   => self%data_c)

      call check(ia >= 1 .and. ia+ma-1 <= m,'set_value_MatLAP_c','row indices out of range')
      call check(ja >= 1 .and. ja+na-1 <= n,'set_value_MatLAP_c','column indices out of range')

      if(type == 'dc') then
         
         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               if(incl) then
                  ac(i,j) = ac(i,j) + a(i-ia+1,j-ja+1)
               else
                  ac(i,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do

      else if(type == 'dsc' .or. type == 'dhp') then

         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               if(i > j) cycle
               if(incl) then
                  ac(i,j) = ac(i,j) + a(i-ia+1,j-ja+1)
               else
                  ac(i,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do

      else if(type == 'bc') then

         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               k = kl + ku + 1 + i - j
               if(k < 1 .or. k > ld) cycle
               if(incl) then
                  ac(k,j) = ac(k,j) + a(i-ia+1,j-ja+1)
               else
                  ac(k,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do         

      else if(type == 'bhp') then
         
         do j = ja,ja+na-1
            do i = ia,ia+ma-1
               k = kd + 1 + i - j
               if(k < 1 .or. k > ld) cycle
               if(incl) then
                  ac(k,j) = ac(k,j) + a(i-ia+1,j-ja+1)
               else
                  ac(k,j) = a(i-ia+1,j-ja+1)                  
               end if
            end do
         end do
         
      else

         call error('set_value_MatLAP_c','unknown matrix type')
         
      end if
      
    end associate
         
    return
  end subroutine set_value_MatLAP_c


  subroutine set_value_MatLAP_single_c(self,i,j,a,inc)
    class(MatLAP), intent(inout) :: self
    integer(i4b), intent(in) :: i,j
    complex(dpc), intent(in) :: a
    logical, intent(in), optional :: inc
    complex(dpc), dimension(1,1) :: al
    al = a
    call self%set(i,1,j,1,al,inc)
    return
  end subroutine set_value_MatLAP_single_c


  

  subroutine factor_MatLAP(self)
    class(MatLAP), intent(inout) :: self
    
    logical :: test
    integer(i4b) :: p,info,lwork
    real(dp), dimension(1) :: work_r_tmp
    real(dp), dimension(:), allocatable :: work_r
    complex(dpc), dimension(1) :: work_c_tmp
    complex(dpc), dimension(:), allocatable :: work_c

    call check(self%allocated,'factor_MatLAP','matrix not allocated')
    
    if(self%factored) return
    
    test =      self%type == 'dr'  &
           .or. self%type == 'dc'  &
           .or. self%type == 'dsr' &
           .or. self%type == 'dsc' &
           .or. self%type == 'dh'  &
           .or. self%type == 'br'  &
           .or. self%type == 'bc' 
    if(test) then
       p = min(self%m,self%n)
       allocate(self%ipiv(p))
    end if
    
    associate(type => self%type,    &
              m    => self%m,       &
              n    => self%n,       &
              kl   => self%kl,      &
              ku   => self%ku,      &
              kd   => self%kd,      &
              ld   => self%ld,      & 
              ar   => self%data_r,  &
              ac   => self%data_c,  &
              ipiv => self%ipiv)
      
      if(type == 'dr') then
         
         call dgetrf(m,n,ar,m,ipiv,info)
         
      else if(type == 'dc') then
         
         call zgetrf(m,n,ac,m,ipiv,info)
         
      else if(type == 'dsr') then
         
         call dsytrf('U',n,ar,n,ipiv,work_r_tmp,-1,info)
         lwork = work_r_tmp(1)
         allocate(work_r(lwork))
         call dsytrf('U',n,self%data_r,n,ipiv,work_r,lwork,info)
         
      else if(type == 'dsc') then
         
         call zsytrf('U',n,ac,n,ipiv,work_c_tmp,-1,info)
         lwork = work_c_tmp(1)
         allocate(work_c(lwork))
         call zsytrf('U',n,ac,n,ipiv,work_c,lwork,info)
         
      else if(type == 'dh') then
         
         call zhetrf('U',n,ac,n,ipiv,work_c_tmp,-1,info)
         lwork = work_c_tmp(1)
         allocate(work_c(lwork))
         call zhetrf('U',n,ac,n,ipiv,work_c,lwork,info)
         
      else if(type == 'dsp') then
         
         call dpotrf('U',n,ar,n,info)
         
      else if(type == 'dhp') then
         
         call zpotrf('U',n,ac,n,info)
         
      else if(type == 'br') then
         
         call dgbtrf(m,n,kl,ku,ar,ld,ipiv,info)
         
      else if(type == 'bc') then
         
         call zgbtrf(m,n,kl,ku,ac,ld,ipiv,info)
         
      else if(type == 'bsp') then
         
         call dpbtrf('U',n,kd,ar,ld,info)
         
      else if(type == 'bhp') then
         
         call zpbtrf('U',n,kd,ac,ld,info)
         
      else
         
         call error('allocate_MatLAP','undefined type')
         
      end if
      
      call check(info == 0,'allocate_MatLAP','problem with factorisation')
      self%factored = .true.
      
    end associate
    return
  end subroutine factor_MatLAP
  
  subroutine backsub_MatLAP_r_single(self,b,trans)
    class(MatLAP), intent(in) :: self
    real(dp), dimension(:), intent(in), target :: b
    character(len=1), intent(in), optional  :: trans
    
    character(len=1) :: transl
    character(len=1), dimension(2) :: trans_op = (/'N','T'/)
    integer(i4b), parameter :: nrhs = 1
    integer(i4b) :: info
    real(dp), dimension(:,:), pointer :: bl

    call check(self%allocated,'backsub_MatLAP_r_single','matrix not allocated')
    call check(self%factored,'backsub_MatLAP_r_single','matrix must be factored')
    
    if(present(trans)) then
       call check(any(trans_op == trans),'backsub_MatLAP_r_single','invalid option')
       transl = trans
    else
       transl = 'N'
    end if
    
    associate(type => self%type,    &
              m    => self%m,       &
              n    => self%n,       &
              kl   => self%kl,      &
              ku   => self%ku,      &
              kd   => self%kd,      &
              ld   => self%ld,      & 
              a    => self%data_r,  &
              ipiv => self%ipiv)

      call check(m == n,'backsub_MatLAP_r_single','matrix must be square')
      call check(n == size(b),'backsub_MatLAP_r_single','b has wrong dimensions')
      
      bl(1:n,1:1) => b
      
      if(type == 'dr') then

         call dgetrs(transl,n,nrhs,a,n,ipiv,bl,n,info)

      else if(type == 'dsr') then

         call dsytrs('U',n,nrhs,a,n,ipiv,bl,n,info)

      else if(type == 'dsp') then

         call dpotrs('U',n,nrhs,a,n,bl,n,info)

      else if(type == 'br') then

         call dgbtrs(transl,n,kl,ku,nrhs,a,ld,ipiv,bl,n,info)

      else if(type == 'bsp') then

         call dpbtrs('U',n,kd,nrhs,a,ld,bl,n,info)
         
      else

         call error('backsub_MatLAP_r_single','unknown matrix type')
         
      end if
      
      call check(info == 0,'backsub_MatLAP_r_single','problem with back substitution')      
      nullify(bl)
      
    end associate
    
    return
  end subroutine backsub_MatLAP_r_single

  
  subroutine backsub_MatLAP_r_many(self,b,trans)
    class(MatLAP), intent(in) :: self
    real(dp), dimension(:,:), intent(in) :: b
    character(len=1), intent(in), optional  :: trans
    
    character(len=1) :: transl
    character(len=1), dimension(2) :: trans_op = (/'N','T'/)
    integer(i4b) :: nrhs,info

    call check(self%allocated,'backsub_MatLAP_r_many','matrix not allocated')
    call check(self%factored,'backsub_MatLAP_r_many','matrix must be factored')
    
    if(present(trans)) then
       call check(any(trans_op == trans),'backsub_MatLAP_r_many','invalid option')
       transl = trans
    else
       transl = 'N'
    end if
    
    associate(type => self%type,    &
              m    => self%m,       &
              n    => self%n,       &
              kl   => self%kl,      &
              ku   => self%ku,      &
              kd   => self%kd,      &
              ld   => self%ld,      & 
              a    => self%data_r,  &
              ipiv => self%ipiv)

      call check(m == n,'backsub_MatLAP_r_many','matrix must be square')
      call check(n == size(b),'backsub_MatLAP_r_many','b has wrong dimensions')
      
      nrhs = size(b,2)
      
      if(type == 'dr') then

         call dgetrs(transl,n,nrhs,a,n,ipiv,b,n,info)

      else if(type == 'dsr') then

         call dsytrs('U',n,nrhs,a,n,ipiv,b,n,info)

      else if(type == 'dsp') then

         call dpotrs('U',n,nrhs,a,n,b,n,info)

      else if(type == 'br') then

         call dgbtrs(transl,n,kl,ku,nrhs,a,ld,ipiv,b,n,info)

      else if(type == 'bsp') then

         call dpbtrs('U',n,kd,nrhs,a,ld,b,n,info)
         
      else

         call error('backsub_MatLAP_r_many','unknown matrix type')
         
      end if
      
      call check(info == 0,'backsub_MatLAP_r_many','problem with back substitution')      
      
    end associate
    
    return
  end subroutine backsub_MatLAP_r_many


  
  subroutine backsub_MatLAP_c_single(self,b,trans)
    class(MatLAP), intent(in) :: self
    complex(dpc), dimension(:), intent(in), target :: b
    character(len=1), intent(in), optional  :: trans
    
    character(len=1) :: transl
    character(len=1), dimension(3) :: trans_op = (/'N','T','C'/)
    integer(i4b), parameter :: nrhs = 1
    integer(i4b) :: info
    complex(dpc), dimension(:,:), pointer :: bl

    call check(self%allocated,'backsub_MatLAP_c_single','matrix not allocated')
    call check(self%factored,'backsub_MatLAP_c_single','matrix must be factored')
    
    if(present(trans)) then
       call check(any(trans_op == trans),'backsub_MatLAP_c_single','invalid option')
       transl = trans
    else
       transl = 'N'
    end if
    
    associate(type => self%type,    &
              m    => self%m,       &
              n    => self%n,       &
              kl   => self%kl,      &
              ku   => self%ku,      &
              kd   => self%kd,      &
              ld   => self%ld,      & 
              a    => self%data_c,  &
              ipiv => self%ipiv)

      call check(m == n,'backsub_MatLAP_c_single','matrix must be square')
      call check(n == size(b),'backsub_MatLAP_c_single','b has wrong dimensions')
      
      bl(1:n,1:1) => b
      
      if(type == 'dc') then

         call zgetrs(transl,n,nrhs,a,n,ipiv,bl,n,info)

      else if(type == 'dsc') then

         call zsytrs('U',n,nrhs,a,n,ipiv,bl,n,info)

      else if(type == 'dhp') then

         call zpotrs('U',n,nrhs,a,n,bl,n,info)

      else if(type == 'bc') then

         call zgbtrs(transl,n,kl,ku,nrhs,a,ld,ipiv,bl,n,info)

      else if(type == 'bhp') then

         call zpbtrs('U',n,kd,nrhs,a,ld,bl,n,info)
         
      else

         call error('backsub_MatLAP_c_single','unknown matrix type')
         
      end if
      
      call check(info == 0,'backsub_MatLAP_c_single','problem with back substitution')      
      nullify(bl)
      
    end associate
    
    return
  end subroutine backsub_MatLAP_c_single

  
  subroutine backsub_MatLAP_c_many(self,b,trans)
    class(MatLAP), intent(in) :: self
    complex(dpc), dimension(:,:), intent(in) :: b
    character(len=1), intent(in), optional  :: trans
    
    character(len=1) :: transl
    character(len=1), dimension(3) :: trans_op = (/'N','T','C'/)
    integer(i4b) :: nrhs,info

    call check(self%allocated,'backsub_MatLAP_c_many','matrix not allocated')
    call check(self%factored,'backsub_MatLAP_c_many','matrix must be factored')
    
    if(present(trans)) then
       call check(any(trans_op == trans),'backsub_MatLAP_c_many','invalid option')
       transl = trans
    else
       transl = 'N'
    end if
    
    associate(type => self%type,    &
              m    => self%m,       &
              n    => self%n,       &
              kl   => self%kl,      &
              ku   => self%ku,      &
              kd   => self%kd,      &
              ld   => self%ld,      & 
              a    => self%data_c,  &
              ipiv => self%ipiv)

      call check(m == n,'backsub_MatLAP_c_many','matrix must be square')
      call check(n == size(b),'backsub_MatLAP_c_many','b has wrong dimensions')
      
      nrhs = size(b,2)
      
      if(type == 'dc') then

         call zgetrs(transl,n,nrhs,a,n,ipiv,b,n,info)

      else if(type == 'dsc') then

         call zsytrs('U',n,nrhs,a,n,ipiv,b,n,info)

      else if(type == 'dhp') then

         call zpotrs('U',n,nrhs,a,n,b,n,info)

      else if(type == 'bc') then

         call zgbtrs(transl,n,kl,ku,nrhs,a,ld,ipiv,b,n,info)

      else if(type == 'bhp') then

         call zpbtrs('U',n,kd,nrhs,a,ld,b,n,info)
         
      else

         call error('backsub_MatLAP_c_many','unknown matrix type')
         
      end if
      
      call check(info == 0,'backsub_MatLAP_c_many','problem with back substitution')      
      
    end associate
    
    return
  end subroutine backsub_MatLAP_c_many
  

  
end module module_LAPACK
