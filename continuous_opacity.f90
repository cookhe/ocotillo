! Module to calculate hydrogen cross sections as described
!   in Gray's (2022) description of continuous opacity in
!   stellar atmospheres.
module ContinuousOpacity

  use Common

  implicit none
  private
  
  public :: calc_hydrogen_stimulated_emission
  public :: get_electron_thomson_scattering,get_theta
  public :: get_hydrogen_ion_bound_free
  public :: calc_opacity_and_albedo
  
contains
!******************************************
  function get_electron_thomson_scattering(n, nHII) result(e_scatter)
    real, intent(in), dimension(nz) :: n, nHII
    real, dimension(nz) ::  e_scatter
    real :: alpha_e = 0.6648e-24 ! coefficient
    
    e_scatter = alpha_e * (nHII / n)  
  
  endfunction get_electron_thomson_scattering
!******************************************
  function get_theta(T) result(theta)
    real, intent(in), dimension(nz) :: T
    real, dimension(nz) ::  theta
    
    theta = 5040./ T
    
  endfunction get_theta
!******************************************
  function get_hydrogen_ion_bound_free(electron_pressure,theta) result(hm_bf_factor)

    real, intent(in), dimension(nz) :: electron_pressure,theta
    real, dimension(nz) :: hm_bf_factor
    
    hm_bf_factor=4.158e-10 * electron_pressure * theta**(5./2) * 10**(0.754 * theta)

  endfunction get_hydrogen_ion_bound_free
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
  subroutine calc_opacity_and_albedo(opacity, albedo)

    real, dimension(nz) :: opacity, albedo

    intent(out)  :: opacity, albedo

    !kappa_H_bf  = specb.xsec_bfHI(np.array([wave_angstrom]), T, m=6)
    !kappa_H_bf *= stim_factor * ionization_factor
    !kappa_H_ff   = specb.xsec_ffH(np.array(wave_angstrom), T)
    !kappa_H_ff  *= stim_factor * ionization_factor
    !!
    !alpha_Hmbf          = specb.xsec_bfHminus(np.array([wave_angstrom]))
    !kappa_Hm_bf  = hm_bf_factor * alpha_Hmbf
    !kappa_Hm_bf *= stim_factor * ionization_factor
    !
    !kappa_Hm_ff  = specb.xsec_ffHminus(ne, np.array([wave_angstrom]), T)
    !kappa_Hm_ff *= ionization_factor

    !kappa_rad = kappa_H_bf + kappa_H_ff + kappa_Hm_bf + kappa_Hm_ff
    !kappa_rad += e_scatter

    !! placeholder until albedo properly calculated
    !albedo = (e_scatter) / (kappa_rad)
    !kappa_rad /= mp
    
    !! placeholder until absorption coefficient properly calculated
    !opacity = kappa_rad * rho

    albedo=1
    opacity=1

    
  endsubroutine calc_opacity_and_albedo
!******************************************
endmodule ContinuousOpacity
