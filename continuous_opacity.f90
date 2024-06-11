! Module to calculate hydrogen cross sections as described
!   in Gray's (2022) description of continuous opacity in
!   stellar atmospheres.
module ContinuousOpacity

  use Common

  implicit none
  private
  
  public :: calc_hydrogen_stimulated_emission,calc_electron_thomson_scattering

contains
!******************************************
subroutine calc_electron_thomson_scattering(n, nHII, e_scatter)
  real, intent(in), dimension(nz) :: n, nHII
  real, intent(out), dimension(nz) ::  e_scatter
  real :: alpha_e = 0.6648e-24 ! coefficient
  
  e_scatter = alpha_e * (nHII / n)  
  
endsubroutine calc_electron_thomson_scattering
!******************************************
subroutine calc_hydrogen_stimulated_emission(wave, theta, stim_factor)
  real, intent(in) :: wave  ! MUST be in angstroms
  real, intent(in), dimension(nz) :: theta  ! theta = 5040./T[K]
  real, intent(out), dimension(nz) :: stim_factor
  real :: chi

  chi = 1.2398e4 / wave
  stim_factor = 1 - 10**(-chi*theta)

endsubroutine calc_hydrogen_stimulated_emission
!******************************************



endmodule ContinuousOpacity