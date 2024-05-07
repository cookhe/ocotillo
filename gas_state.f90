! Module to calculate the gas state
module gas_state

  use common


  implicit none
  private

  public :: hydrogen_ion_frac

  real :: fully_ionized_T

  namelist /gas_state_input/ fully_ionized_T

contains
!******************************************
subroutine hydrogen_ion_frac(rho,T,NHII_NHINHII)
  real, dimension(nz), intent(in) :: rho,T
  real, dimension(nz), intent(out) :: NHII_NHINHII
  real, dimension(nz) :: constants,niine_ni,C,index
  
  ! calculate the Saha equation (relative fraction of adjacent ions)
  constants = ((sqrt(2*pi*me*k_cgs*T)/h_planck)**3)
  index = hydrogen_ionization_eV/(k_eV*T)
  ! index = 100.
  ! print*, 'E/k_eV*T = ', index
  niine_ni =  constants * exp(-index)
    
  C = mp / rho * niine_ni

  NHII_NHINHII = (sqrt(C**2 + 4*C) - C)/2.

  ! Limits due to machine precision mean 4C becomes unresolved next to C**2,
  ! so manually set the ionization fraction to 1 for temperatures above 20000.
  ! This is appropriate for our domain where the density is always less than
  ! ~2e-9 g cm^-3. Would not be appropriate for higher densities.
  ! print*, where T > fully_ionized_T
  NHII_NHINHII = merge(1.0, NHII_NHINHII, T>fully_ionized_T)
  
end subroutine hydrogen_ion_frac

end module gas_state
