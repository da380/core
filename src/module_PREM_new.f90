module module_PREM_new

  use module_constants
  use module_physical_constants
  use module_util
  use module_spherical_model_new
  implicit none


  type, extends(density_model) :: PREM_density_model
     logical :: ocean
   contains
     procedure :: deallocate => deallocate_PREM_density_model
     procedure :: rho  => rho_PREM_density_model
     procedure :: drho => drho_PREM_density_model
     procedure :: phi  => phi_PREM_density_model
     procedure :: g    => g_PREM_density_model
     procedure :: ep   => ep_PREM_density_model
     procedure :: mass => mass_PREM_density_model
     procedure :: gs => mass_PREM_density_model
     procedure :: I1 => mass_PREM_density_model
     procedure :: I2 => mass_PREM_density_model
     procedure :: I3 => mass_PREM_density_model          
  end type PREM_density_model

  interface PREM_density_model
     procedure :: build_PREM_density_model
  end interface PREM_density_model
  

  real(dp), dimension(2,13), parameter :: r_PREM = reshape((/   0.0e3_dp, 1221.5e3_dp, &
                                                             1221.5e3_dp, 3480.0e3_dp, & 
                                                             3480.0e3_dp, 3630.0e3_dp, &
                                                             3630.0e3_dp, 5600.0e3_dp, &
                                                             5600.0e3_dp, 5701.0e3_dp, &
                                                             5701.0e3_dp, 5771.0e3_dp, &
                                                             5771.0e3_dp, 5971.0e3_dp, & 
                                                             5971.0e3_dp, 6151.0e3_dp, & 
                                                             6151.0e3_dp, 6291.0e3_dp, &
                                                             6291.0e3_dp, 6346.6e3_dp, &
                                                             6346.6e3_dp, 6356.0e3_dp, &
                                                             6356.0e3_dp, 6368.0e3_dp, &
                                                             6368.0e3_dp, 6371.0e3_dp/)/length_norm,(/2,13/))

                                                             
  real(dp), dimension(4,13), parameter :: rho_coef_PREM = reshape((/ 13.08850e3_dp, 0.0e3_dp,-8.83810e3_dp, 0.0e3_dp,            &
                                                                12.58150e3_dp, -1.26380e3_dp, -3.64260e3_dp, -5.52810e3_dp, &
                                                                7.95650e3_dp, -6.47610e3_dp,  5.52830e3_dp, -3.08070e3_dp,  &
                                                                7.95650e3_dp, -6.47610e3_dp,  5.52830e3_dp, -3.08070e3_dp,  &
                                                                7.95650e3_dp, -6.47610e3_dp,  5.52830e3_dp, -3.08070e3_dp,  &
                                                                5.31970e3_dp, -1.48360e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                11.24940e3_dp, -8.02980e3_dp , 0.0e3_dp, 0.0e3_dp,          &
                                                                7.10890e3_dp, -3.80450e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                2.69100e3_dp,  0.69240e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                2.69100e3_dp,  0.69240e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                2.90000e3_dp, 0.0e3_dp, 0.0e3_dp, 0.0e3_dp,                 &
                                                                2.60000e3_dp, 0.0e3_dp, 0.0e3_dp, 0.0e3_dp,                 &
                                                                1.02000e3_dp, 0.0e3_dp, 0.0e3_dp, 0.0e3_dp/)/density_norm,  &
                                                                (/4,13/))

  real(dp), dimension(4,13), parameter :: vpv_coef_PREM = reshape((/ 11.26220e3_dp,  0.0e3_dp, -6.36400e3_dp, 0.0e3_dp,          &
                                                                11.04870e3_dp, -4.03620e3_dp,  4.8023e3_dp, -13.57320e3_dp, &
                                                                15.38910e3_dp, -5.31810e3_dp,  5.52420e3_dp, -2.55140e3_dp, &
                                                                24.9520e3_dp, -40.4673e3_dp,  51.4832e3_dp, -26.64190e3_dp, &
                                                                29.2766e3_dp, -23.6027e3_dp,   5.52420e3_dp, -2.55140e3_dp, &
                                                                19.09570e3_dp, -9.86720e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                39.7027e3_dp, -32.61660e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                20.3926e3_dp, -12.25690e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                0.83170e3_dp,  7.21800e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                0.83170e3_dp,  7.21800e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                6.80000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp,                &
                                                                5.80000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp,                &
                                                                1.45000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp/)/velocity_norm,&
                                                                (/4,13/))
  
  
  real(dp), dimension(4,13), parameter :: vsv_coef_PREM = reshape((/ 3.66780e3_dp,  0.0e3_dp, -4.44750e3_dp, 0.0e3_dp,           &
                                                                0.0e3_dp,0.0e3_dp,0.0e3_dp,0.0e3_dp,                        &
                                                                6.92540e3_dp,  1.46720e3_dp, -2.08340e3_dp,  0.97830e3_dp,  &
                                                                11.1671e3_dp, -13.7818e3_dp,  17.4575e3_dp,  -9.2777e3_dp,  & 
                                                                22.3459e3_dp, -17.2473e3_dp,  -2.08340e3_dp,  0.97830e3_dp, & 
                                                                9.98390e3_dp, -4.93240e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                22.3512e3_dp, -18.58560e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                8.94960e3_dp, -4.45970e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                5.85820e3_dp, -1.46780e3_dp , 0.0e3_dp, 0.0e3_dp,           &
                                                                5.85820e3_dp, -1.46780e3_dp , 0.0e3_dp, 0.0e3_dp,           &
                                                                3.90000e3_dp, 0.0e3_dp , 0.0e3_dp, 0.0e3_dp,                &
                                                                3.20000e3_dp, 0.0e3_dp , 0.0e3_dp, 0.0e3_dp,                &
                                                                0.0e3_dp,0.0e3_dp,0.0e3_dp,0.0e3_dp/)/velocity_norm,(/4,13/))

  real(dp), dimension(4,13), parameter :: qk_coef_PREM  = reshape((/ 1327.7_dp, 0.0_dp, 0.0_dp, 0.0_dp,                  &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/),(/4,13/))



  real(dp), dimension(4,13), parameter :: qm_coef_PREM  = reshape((/ 84.6_dp, 0.0_dp, 0.0_dp, 0.0_dp,                    &
                                                                0.0_dp,  0.0_dp, 0.0_dp, 0.0_dp,                    &
                                                                312.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                312.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                312.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                143.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                143.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                143.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                80.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                    &
                                                                600.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                600.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                600.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/),(/4,13/))


  real(dp), dimension(4,13), parameter :: vph_coef_PREM = reshape((/ 11.26220e3_dp,  0.0e3_dp, -6.36400e3_dp, 0.0e3_dp,          &
                                                                11.04870e3_dp, -4.03620e3_dp,  4.8023e3_dp, -13.57320e3_dp, &
                                                                15.38910e3_dp, -5.31810e3_dp,  5.52420e3_dp, -2.55140e3_dp, &
                                                                24.9520e3_dp, -40.4673e3_dp,  51.4832e3_dp, -26.64190e3_dp, &
                                                                29.2766e3_dp, -23.6027e3_dp,   5.52420e3_dp, -2.55140e3_dp, &
                                                                19.09570e3_dp, -9.86720e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                39.7027e3_dp, -32.61660e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                20.3926e3_dp, -12.25690e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                3.59080e3_dp,  4.61720e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                3.59080e3_dp,  4.61720e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                6.80000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp,                &
                                                                5.80000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp,                &
                                                                1.45000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp/)/velocity_norm,&
                                                                (/4,13/))

  real(dp), dimension(4,13), parameter :: vsh_coef_PREM = reshape((/ 3.66780e3_dp,  0.0e3_dp, -4.44750e3_dp, 0.0e3_dp,           &
                                                                0.0e3_dp,0.0e3_dp,0.0e3_dp,0.0e3_dp,                        &
                                                                6.92540e3_dp,  1.46720e3_dp, -2.08340e3_dp,  0.97830e3_dp,  &
                                                                11.1671e3_dp, -13.7818e3_dp,  17.4575e3_dp,  -9.2777e3_dp,  & 
                                                                22.3459e3_dp, -17.2473e3_dp,  -2.08340e3_dp,  0.97830e3_dp, & 
                                                                9.98390e3_dp, -4.93240e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                22.3512e3_dp, -18.58560e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                8.94960e3_dp, -4.45970e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                -1.08390e3_dp,  5.71760e3_dp , 0.0e3_dp, 0.0e3_dp,          &
                                                                -1.08390e3_dp,  5.71760e3_dp , 0.0e3_dp, 0.0e3_dp,          &
                                                                3.90000e3_dp, 0.0e3_dp , 0.0e3_dp, 0.0e3_dp,                &
                                                                3.20000e3_dp, 0.0e3_dp , 0.0e3_dp, 0.0e3_dp,                &
                                                                0.0e3_dp,0.0e3_dp,0.0e3_dp,0.0e3_dp/)/velocity_norm,(/4,13/))


  real(dp), dimension(4,13), parameter :: eta_coef_PREM = reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                3.36870_dp, -2.47780_dp ,0.0_dp,0.0_dp,             &
                                                                3.36870_dp, -2.47780_dp ,0.0_dp,0.0_dp,             &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/),(/4,13/))



contains

  !=============================================================================!
  !                        polynomial evaluation functions                      !
  !=============================================================================!

  
  real(dp) function rho_PREM(i,r) result(rho)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    rho = poly(4,rho_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function rho_PREM


  real(dp) function drho_PREM(i,r) result(drho)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    drho = dpoly(4,rho_coef_PREM(:,i),r/r_PREM(2,13))/r_PREM(2,13)
    return
  end function drho_PREM

  
  real(dp) function  vpv_PREM(i,r) result(vpv)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    vpv = poly(4,vpv_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function vpv_PREM

  
  real(dp) function  vph_PREM(i,r) result(vph)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    vph = poly(4,vph_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function vph_PREM


  real(dp) function  vsv_PREM(i,r) result(vsv)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    vsv = poly(4,vsv_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function vsv_PREM

  
  real(dp) function  vsh_PREM(i,r) result(vsh)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    vsh = poly(4,vsh_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function vsh_PREM


  real(dp) function  qk_PREM(i,r) result(qk)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    qk = poly(4,qk_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function qk_PREM


  real(dp) function  qm_PREM(i,r) result(qm)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    qm = poly(4,qm_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function qm_PREM


  real(dp) function  eta_PREM(i,r) result(eta)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    eta = poly(4,eta_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function eta_PREM


  !=============================================================================!
  !                            density model procedures                         !
  !=============================================================================!


  type(PREM_density_model) function build_PREM_density_model(ocean) result(model)
    logical, intent(in), optional :: ocean

    integer(i4b) :: i,n
    
    if(present(ocean)) model%ocean = ocean
    if(model%ocean) then
       n = 13
    else
       n = 12
    end if
    model%n = n
    
    allocate(model%r(n+1))
    do i = 1,n
       model%r(i) = r_PREM(1,i)
    end do
    model%r(n+1) = r_PREM(2,n)
    
    return
  end function build_PREM_density_model
  

  subroutine deallocate_PREM_density_model(self)
    class(PREM_density_model), intent(inout) :: self
    deallocate(self%r)
    return
  end subroutine deallocate_PREM_density_model

  real(dp) function rho_PREM_density_model(self,i,r) result(rho)
    class(PREM_density_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    rho = rho_PREM(i,r)
    return
  end function rho_PREM_density_model


  real(dp) function drho_PREM_density_model(self,i,r) result(drho)
    class(PREM_density_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    drho = drho_PREM(i,r)
    return
  end function drho_PREM_density_model


  real(dp) function phi_PREM_density_model(self,i,r) result(phi)
    class(PREM_density_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    phi = 0.0_dp
    return
  end function phi_PREM_density_model


  real(dp) function g_PREM_density_model(self,i,r) result(g)
    class(PREM_density_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    g = 0.0_dp
    return
  end function g_PREM_density_model


  real(dp) function ep_PREM_density_model(self,i,r) result(ep)
    class(PREM_density_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    ep = 0.0_dp
    return
  end function ep_PREM_density_model


  real(dp) function mass_PREM_density_model(self) result(mass)
    class(PREM_density_model), intent(in) :: self
    mass = 0.0_dp
    return
  end function mass_PREM_density_model


  real(dp) function gs_PREM_density_model(self) result(gs)
    class(PREM_density_model), intent(in) :: self
    gs = 0.0_dp
    return
  end function gs_PREM_density_model


  real(dp) function I1_PREM_density_model(self) result(I1)
    class(PREM_density_model), intent(in) :: self
    I1 = 0.0_dp
    return
  end function I1_PREM_density_model


  real(dp) function I2_PREM_density_model(self) result(I2)
    class(PREM_density_model), intent(in) :: self
    I2 = 0.0_dp
    return
  end function I2_PREM_density_model


  real(dp) function I3_PREM_density_model(self) result(I3)
    class(PREM_density_model), intent(in) :: self
    I3 = 0.0_dp
    return
  end function I3_PREM_density_model
  
  
end module module_PREM_new
