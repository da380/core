module module_spherical_model_new

  use module_constants
  use module_error
  use module_util
  implicit none


  !======================================================!
  !        set the base class for spherical models       !
  !======================================================!
  
  type, abstract :: density_model
     integer(i4b) :: n
     real(dp), dimension(:), allocatable :: r
   contains
     procedure(deallocate_model), deferred :: deallocate 
     procedure :: r1 => r1_spherical_model
     procedure :: r2 => r2_spherical_model
     procedure :: find_layer => find_layer_spherical_model
     procedure(rho_function), deferred :: rho
     procedure(rho_function), deferred :: drho
     procedure(rho_function), deferred :: phi
     procedure(rho_function), deferred :: g
     procedure(rho_function), deferred :: ep     
     procedure(mass_function), deferred :: mass
     procedure(mass_function), deferred :: gs
     procedure(mass_function), deferred :: I1
     procedure(mass_function), deferred :: I2
     procedure(mass_function), deferred :: I3     
  end type density_model

  abstract interface

     subroutine deallocate_model(self)
       import density_model
       class(density_model), intent(inout) :: self
     end subroutine deallocate_model
     
     real(dp) function rho_function(self,i,r) 
       use module_constants
       import density_model
       class(density_model), intent(in) :: self
       integer(i4b), intent(in) :: i
       real(dp), intent(in) :: r
     end function rho_function

     real(dp) function mass_function(self) 
       use module_constants
       import density_model
       class(density_model), intent(in) :: self       
     end function mass_function
     
  end interface
  
contains

  real(dp) function r1_spherical_model(self) result(r1)
    class(density_model), intent(in) :: self
    r1 = self%r(1)
    return
  end function r1_spherical_model

  real(dp) function r2_spherical_model(self) result(r2)
    class(density_model), intent(in) :: self
    r2 = self%r(self%n+1)
    return
  end function r2_spherical_model

  integer(i4b) function find_layer_spherical_model(self,r) result(i)
    class(density_model), intent(in) :: self
    real(dp), intent(in) :: r
    i = bisect_list(self%r,r)
    return
  end function find_layer_spherical_model
  
  
end module module_spherical_model_new
