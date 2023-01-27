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
  !             hermitian matrix declarations            !
  !------------------------------------------------------!

  interface allocate_matrix_h
     procedure :: allocate_matrix_s_real
     procedure :: allocate_matrix_s_complex
  end interface allocate_matrix_h

  interface factorise_matrix_h
     procedure :: factorise_matrix_s_real
     procedure :: factorise_matrix_h_complex     
  end interface factorise_matrix_h
  
  interface backsub_matrix_h
     procedure :: backsub_matrix_s_real_single
     procedure :: backsub_matrix_h_complex_many     
  end interface backsub_matrix_h


  interface eigenvalue_matrix_h
  end interface eigenvalue_matrix_h
  
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
  
  interface allocate_matrix_bhp
     procedure :: allocate_matrix_bhp_real
     procedure :: allocate_matrix_bhp_complex
  end interface allocate_matrix_bhp
  
  interface factorise_matrix_bhp     
     procedure :: factorise_matrix_bhp_real
     procedure :: factorise_matrix_bhp_complex
  end interface factorise_matrix_bhp

  interface backsub_matrix_bhp
     procedure :: backsub_matrix_bhp_real_single
     procedure :: backsub_matrix_bhp_real_many
     procedure :: backsub_matrix_bhp_complex_single
     procedure :: backsub_matrix_bhp_complex_many
  end interface backsub_matrix_bhp


  !-----------------------------------------------!
  !               other declarations              !
  !-----------------------------------------------!

  interface matrix_s_to_matrix
     procedure ::  matrix_s_to_matrix_real
     procedure ::  matrix_s_to_matrix_complex
  end interface matrix_s_to_matrix

  interface matrix_h_to_matrix
     procedure ::  matrix_s_to_matrix_real
     procedure ::  matrix_h_to_matrix_complex
  end interface matrix_h_to_matrix

  interface matrix_b_to_matrix
     procedure ::  matrix_b_to_matrix_real
     procedure ::  matrix_b_to_matrix_complex
  end interface matrix_b_to_matrix

  interface matrix_bhp_to_matrix
     procedure ::  matrix_bhp_to_matrix_real
     procedure ::  matrix_bhp_to_matrix_complex
  end interface matrix_bhp_to_matrix

  interface matrix_bhp_to_matrix_s
     procedure ::  matrix_bhp_to_matrix_s_real
     procedure ::  matrix_bhp_to_matrix_s_complex
  end interface matrix_bhp_to_matrix_s
  
  interface matrix_bhp_to_matrix_b
     procedure ::  matrix_bhp_to_matrix_b_real
     procedure ::  matrix_bhp_to_matrix_b_complex
  end interface matrix_bhp_to_matrix_b
  
contains


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
    call check(info == 0,'backsub_matrix_s_real_single', &
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
    call check(info == 0,'backsub_matrix_s_real_many', &
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
    call check(info == 0,'backsub_matrix_s_complex_single', &
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
    call check(info == 0,'backsub_matrix_s_complex_many', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_s_complex_many


  !--------------------------------------------------------------------!
  !                  routines for hermitian matrices                   !
  !--------------------------------------------------------------------!
  
  
  subroutine factorise_matrix_h_complex(n,a,ipiv)
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
    call zhetrf('U',n,a,n,ipiv,work,lwork,info)
    call check(info == 0,'factorise_matrix_h_complex','problem with factorisation')
    return
  end subroutine factorise_matrix_h_complex


  subroutine backsub_matrix_h_complex_single(n,a,ipiv,b)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    complex(dpc), dimension(n), intent(inout), target :: b
    integer(i4b) :: info
    complex(dpc), dimension(:,:), pointer :: bl
    bl(1:n,1:1) => b
    call zhetrs('U',n,1,a,n,ipiv,bl,n,info) 
    call check(info == 0,'backsub_matrix_h_complex_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_h_complex_single

 
  subroutine backsub_matrix_h_complex_many(n,a,ipiv,nrhs,b)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n,n), intent(inout) :: a
    integer(i4b), dimension(n), intent(in) :: ipiv
    integer(i4b), intent(in) :: nrhs
    complex(dpc), dimension(n,nrhs), intent(inout) :: b
    integer(i4b) :: info
    call zhetrs('U',n,nrhs,a,n,ipiv,b,n,info) 
    call check(info == 0,'backsub_matrix_h_complex_many', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_h_complex_many

  
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
  

  subroutine allocate_matrix_bhp_real(n,kd,a)
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
  end subroutine allocate_matrix_bhp_real

  subroutine allocate_matrix_bhp_complex(n,kd,a)
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
  end subroutine allocate_matrix_bhp_complex

  integer(i4b) pure function row_index_matrix_bhp(kd,i,j) result(k)
    integer(i4b), intent(in) :: kd,i,j
    k = kd+1+i-j
    return
  end function row_index_matrix_bhp


  subroutine factorise_matrix_bhp_real(n,kd,a)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(kd+1,n), intent(inout) :: a
    integer(i4b) :: ld,info
    ld = kd+1
    call dpbtrf	('U',n,kd,a,ld,info)
    call check(info == 0,'factorise_matrix_bhp_real','problem with decompsition')
    return
  end subroutine factorise_matrix_bhp_real


  subroutine factorise_matrix_bhp_complex(n,kd,a)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(kd+1,n), intent(inout) :: a
    integer(i4b) :: ld,info
    ld = kd+1
    call zpbtrf	('U',n,kd,a,ld,info)
    call check(info == 0,'factorise_matrix_bhp_complex','problem with decompsition')
    return
  end subroutine factorise_matrix_bhp_complex

  
  subroutine backsub_matrix_bhp_real_single(n,kd,a,b)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(kd+1,n), intent(inout) :: a
    real(dp), dimension(n), intent(inout), target :: b
    integer(i4b) :: ld,info
    real(dp), dimension(:,:), pointer :: bl
    ld = kd+1
    bl(1:n,1:1) => b
    call dpbtrs('U',n,kd,1,a,ld,bl,n,info)
    call check(info == 0,'backsub_matrix_bhp_real_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_bhp_real_single


  subroutine backsub_matrix_bhp_real_many(n,kd,a,nrhs,b)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(kd+1,n), intent(inout) :: a
    integer(i4b), intent(in) :: nrhs
    real(dp), dimension(n,nrhs), intent(inout) :: b
    integer(i4b) :: ld,info
    ld = kd+1
    call dpbtrs('U',n,kd,nrhs,a,ld,b,n,info)
    call check(info == 0,'backsub_matrix_bhp_real_many', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_bhp_real_many


  subroutine backsub_matrix_bhp_complex_single(n,kd,a,b)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(kd+1,n), intent(inout) :: a
    complex(dpc), dimension(n), intent(inout), target :: b
    integer(i4b) :: ld,info
    complex(dpc), dimension(:,:), pointer :: bl
    ld = kd+1
    bl(1:n,1:1) => b
    call zpbtrs('U',n,kd,1,a,ld,bl,n,info)
    call check(info == 0,'backsub_matrix_bhp_complex_single', &
                         'problem with backsubtitution')
    nullify(bl)
    return
  end subroutine backsub_matrix_bhp_complex_single


  subroutine backsub_matrix_bhp_complex_many(n,kd,a,nrhs,b)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(kd+1,n), intent(inout) :: a
    integer(i4b), intent(in) :: nrhs
    complex(dpc), dimension(n,nrhs), intent(inout) :: b
    integer(i4b) :: ld,info
    ld = kd+1
    call zpbtrs('U',n,kd,nrhs,a,ld,b,n,info)
    call check(info == 0,'backsub_matrix_bhp_complex_many', &
                         'problem with backsubtitution')
    return
  end subroutine backsub_matrix_bhp_complex_many


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


  subroutine matrix_h_to_matrix_complex(n,a,b)
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n,n), intent(in) :: a
    complex(dpc), dimension(n,n), intent(out) :: b    
    integer(i4b) :: i,j
    do j = 1,n
       do i = 1,j-1
          b(i,j) = a(i,j)
          b(j,i) = conjg(a(i,j))
       end do
       b(j,j) = a(j,j)
       end do
    return
  end subroutine matrix_h_to_matrix_complex


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


  subroutine matrix_bhp_to_matrix_real(n,kd,a,b)
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
  end subroutine matrix_bhp_to_matrix_real


  subroutine matrix_bhp_to_matrix_complex(n,kd,a,b)
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
  end subroutine matrix_bhp_to_matrix_complex


  subroutine matrix_bhp_to_matrix_s_real(n,kd,a,b)
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
  end subroutine matrix_bhp_to_matrix_s_real


  subroutine matrix_bhp_to_matrix_s_complex(n,kd,a,b)
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
  end subroutine matrix_bhp_to_matrix_s_complex
  

  subroutine matrix_bhp_to_matrix_b_real(n,kd,a,kl,ku,b)
    integer(i4b), intent(in) :: n,kd
    real(dp), dimension(kd+1,n), intent(in) :: a
    integer(i4b), intent(in) :: kl,ku
    real(dp), dimension(2*kl+ku+1,n), intent(out) :: b
    integer(i4b) :: i,j,k1,k2    
    call check(kl >= kd,'matrix_bhp_to_matrix_b_real', &
                        'lower band width too small')
    call check(ku >= kd,'matrix_bhp_to_matrix_b_real', &
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
  end subroutine matrix_bhp_to_matrix_b_real


  subroutine matrix_bhp_to_matrix_b_complex(n,kd,a,kl,ku,b)
    integer(i4b), intent(in) :: n,kd
    complex(dpc), dimension(kd+1,n), intent(in) :: a
    integer(i4b), intent(in) :: kl,ku
    complex(dpc), dimension(2*kl+ku+1,n), intent(out) :: b
    integer(i4b) :: i,j,k1,k2    
    call check(kl >= kd,'matrix_bhp_to_matrix_b_complex', &
                        'lower band width too small')
    call check(ku >= kd,'matrix_bhp_to_matrix_b_complex', &
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
  end subroutine matrix_bhp_to_matrix_b_complex




  
end module module_LAPACK
