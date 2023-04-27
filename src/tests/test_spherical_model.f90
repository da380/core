program test_spherical_model

  use module_constants
  use module_physical_constants
  use module_spherical_model_new
  use module_PREM_new

  implicit none
  
  integer(i4b) :: i,j,io,m
  real(dp) :: r1,r2,r,dr,dr0
  class(density_model), allocatable :: model

  model = PREM_density_model(ocean=.true.)
  
  dr0 = 0.01_dp

  open(newunit = io,file='test_spherical_model.out')
  do i = 1,model%n

     r1 = model%r(i)
     r2 = model%r(i+1)
     m = max(int((r2-r1)/dr0),2)
     dr = (r2-r1)/(m-1)
     do j = 1,m
        r = r1+(j-1)*dr
        write(io,*) r*length_norm,model%rho(i,r)*density_norm, &
                    model%drho(i,r)*density_norm/length_norm
     end do
     
  end do
  close(io)
  
  
  
  
end program test_spherical_model


