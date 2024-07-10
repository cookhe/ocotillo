! Module to calculate hydrogen cross sections as described
!   in Gray's (2022) description of continuous opacity in
!   stellar atmospheres.
module ContinuousOpacity

  use Common

  implicit none
  private

  public :: get_hydrogen_stimulated_emission
  public :: get_electron_thomson_scattering
  public :: get_hydrogen_ion_bound_free
  public :: calc_opacity_and_albedo
  public :: grey_parameters
  public :: pre_calc_opacity_quantities
  public :: get_source_function
  public :: calc_wavelength

  real :: switch_ionfraction=1e-2
  real, dimension(7) :: a_coeff
  real, dimension(3,5) :: b_coeff
  real :: AHbf = 1.0449e-26
  real :: Rangstrom=1.0968e-3 ! Rcm = 2 * pi**2 * me * e**4 / (h**3 * c) 
  real :: alpha0=1.0443e-26     ! absorption per electron
  real, dimension(nw,3) :: f_coeff
  integer, parameter :: mdim=6
  real, dimension(mdim) :: n1_array,chi_n
  real, dimension(nw,mdim) :: g
  real :: chi_m
  real, dimension(nw) :: chi,chi1,nu,wa,wa3,gff_factor,wcm1,w1cm5
  real, dimension(nw) :: s_coeff
  real :: damping_factor_const,source_function_const
  
contains
!************************************************************************************
  subroutine calc_wavelength(w1,w0,wa)
!
    real, dimension(nw), intent(out) :: wa
    real, intent(in) :: w1,w0
    real :: dw
    integer :: i
!
    if (nw >1 ) then
      dw = (w1-w0)/(nw-1)
      do i=1,nw
        !wavelengths in angstrom
        wa(i) = w0 + (i-1)*dw
      enddo
    else
      wa=w0
    endif
!     
    print*, 'waves_angstrom min/max', minval(wa),maxval(wa)
    print*, 'number of waves', nw
!
  endsubroutine calc_wavelength
!************************************************************************************
  subroutine pre_calc_opacity_quantities(waves_angstrom)

    real, dimension(nw) :: waves_angstrom,wcm,lgwave
    integer :: n,iw,i,j
!
    damping_factor_const = h_planck*c_light_cgs*k1_cgs
    source_function_const = 2*h_planck*c_light_cgs**2
!    
    wa = waves_angstrom
!    
    wcm  = wa*1d-8 !wavelengths in cm
    wcm1 = 1./wcm
    nu = c_light_cgs*wcm1
!    
    chi  = 1.2398e4/wa
    chi1 = 1./chi
!
    wa3 = wa**3
    w1cm5 = wcm1**5
    gff_factor = 0.3456 * (wa * Rangstrom)**(-one_third)
!
    a_coeff = (/+1.99654,&
         -1.18267e-6,&
         +2.62423e-7,&
         -4.40524e-11,&
         +3.23992e-15,&
         -1.39568e-19,&
         +2.78701e-24/)
    ! loop through each wavelength and perform the summation
    s_coeff=0
    do i=1,7
      s_coeff = s_coeff + 1e-17*a_coeff(i)*wa**(i-1)
    enddo
!
    b_coeff = transpose(reshape(                                   &
         (/-2.276300,-1.685000,+0.766610,-0.053356,+0.000000,&
         +15.28270,-9.284600,+1.993810,-0.142631,+0.000000,  &
         -197.789,+190.266,-67.9775,+10.6913,-0.62515/),     &
         (/ size(b_coeff, 2), size(b_coeff, 1) /)))
!
    do n=1,mdim
      n1_array(n)=1./n
      chi_n(n) = RydbergEnergy * (1. - 1.*n1_array(n)**2) 
      do iw=1,nw
        call gaunt(n, wcm(iw), g(iw,n))
      enddo
    enddo    
    ! excitation energy of state level m
    chi_m = RydbergEnergy * (1. - n1_array(mdim)**2)
!
    lgwave=log10(wa)
    do iw=1,nw
      do i=1,3
        f_coeff(iw,i)=0.
        do j=1,5
          f_coeff(iw,i) = f_coeff(iw,i) + lgwave(iw)**(j-1) * b_coeff(i,j)
        enddo
      enddo
    enddo
!
  endsubroutine pre_calc_opacity_quantities
!************************************************************************************
  function get_electron_thomson_scattering(p) result(e_scatter)

    real, dimension(nz) ::  e_scatter
    real :: alpha_e = 0.6648e-24 ! coefficient
    type (pillar_case) :: p
    
    e_scatter = alpha_e * p%nHII * p%inv_number_density

  endfunction get_electron_thomson_scattering
!************************************************************************************
  function get_hydrogen_ion_bound_free(p) result(hm_bf_factor)

    real, dimension(nz) :: hm_bf_factor
    type (pillar_case) :: p
    
    hm_bf_factor=4.158e-10 * p%electron_pressure * p%theta**2.5 * 10**(0.754 * p%theta)

  endfunction get_hydrogen_ion_bound_free
!************************************************************************************
  function get_hydrogen_stimulated_emission(theta,iw) result(stim_factor)

    real, intent(in), dimension(nz) :: theta  ! theta = 5040./T[K]
    real, dimension(nz) :: stim_factor
    integer :: iw
    
    stim_factor = 1 - 10**(-chi(iw)*theta)

  endfunction get_hydrogen_stimulated_emission
!************************************************************************************
  subroutine calc_opacity_and_albedo(p,iw)

    real, dimension(nz) :: kappa_rad,kappa_H_bf,kappa_Hm_bf,kappa_Hm_ff
    
    integer :: iw
    type (pillar_case) :: p

    p%stim_factor = get_hydrogen_stimulated_emission(p%theta,iw)
    
    kappa_H_bf  = get_kappa_H_bf(p,iw)
    kappa_Hm_bf = get_kappa_Hm_bf(p,iw)
    kappa_Hm_ff = get_kappa_Hm_ff(p,iw)

    kappa_rad = (kappa_H_bf + kappa_Hm_bf + kappa_Hm_ff + p%e_scatter)*mp1
    p%opacity = kappa_rad * p%rho
!
! This subroutine changes p%opacity and kappa_rad
!
    call calc_kappa_H_ff(p,kappa_rad,iw)

    p%albedo = p%e_scatter * mp1 / kappa_rad
    p%opacity1 = 1./p%opacity
    p%opacity2 = p%opacity**2
!
    p%source_function = get_source_function(p%T1,iw)
!    
  endsubroutine calc_opacity_and_albedo
!************************************************************************************
  function get_kappa_H_bf(p, iw) result(kappa_H_bf)
!    
!    """Cross section of bound-free hydrogen. Sums over the first 
!    1 to m-1 excitation states, where m is the principal quantum
!    number and uses the Unsold integral approximation for levels
!    m to infinity. Returns units of cm^2 per neutral hydrogen atom.
!
!    Parameters
!    ----------
!    waves : float
!        Wavelength in Angstroms
!    temp : float
!        Temperature in Kelvin
!    m : int, optional
!        Principal quantum number at which to begin the Unsold
!        integral approximation, by default 6
!
!    Returns
!    -------
!    float
!        Cross section of bound-free hydrogen in units of cm^2 per
!        neutral hydrogen atom
!    """
!    
    real, dimension(nz) :: sm, ktemp, ktemp1, unsold, kappa_H_bf
    integer :: n, iw
    type (pillar_case) :: p
!

    ! sum for the first m-1 excitation states
    sm = 0.
    ktemp = k_cgs*p%T
    ktemp1 = k1_cgs*p%T1
    
    do n=1,mdim
      sm = sm + g(iw,n)*n1_array(n)**3 * exp(-chi_n(n)*ktemp1)
    enddo
    
    ! Unsold approximation integral
    unsold = .5*ktemp*RydbergEnergy1 * ( exp(-chi_m*ktemp1) - exp(-RydbergEnergy*ktemp1) )

    kappa_H_bf = AHbf * p%stim_factor * p%ionization_factor * wa3(iw) * (unsold + sm)
    
  endfunction get_kappa_H_bf
!************************************************************************************
  subroutine calc_kappa_H_ff(p,kappa_rad,iw)
!
!    """Hydrogen free-free absorption coefficient.
!    
!    Units cm^2 per neutral hydrogen atom.
!    """
    real, dimension(nz) :: kappa_rad
    real :: g_ff,energyfactor,kappa_H_ff,opacity_bremsstrahlung
    integer :: i,iw
    type (pillar_case) :: p
!
    intent(inout) :: kappa_rad
!    
    do i=1,nz
       if ((1-p%NHII_NHINHII(i)) .gt. switch_ionfraction) then
          ! Use Gray 2022 function that depends on hydrogen's ionization state.
          g_ff = 1 + gff_factor(iw) * (p%theta1(i)*log10e*chi1(iw) + 0.5)
          energyfactor = .5*log10e*p%theta1(i)*Iev1 * 10**(-p%theta(i)*Iev)
          kappa_H_ff = p%stim_factor(i) * p%ionization_factor(i) * alpha0 * wa3(iw) * g_ff * energyfactor
          kappa_rad(i) = kappa_rad(i) + kappa_H_ff * mp1 ! cm^2 / g
          p%opacity(i) = p%opacity(i) + kappa_rad(i) * p%rho(i) ! 1/cm
       else
          ! Use espression for fully ionized gas.
          !waves_cm = waves / 1e8 ! convert angstroms to centimeters
          !nu_Hz = c*1e8 / waves
          call bremsstrahlung_absorptionCoeff(p,i,nu(iw),1.0,opacity_bremsstrahlung)
          p%opacity(i) = p%opacity(i) + opacity_bremsstrahlung
          kappa_rad(i) = kappa_rad(i) + p%opacity(i) * p%rho1(i) ! cm^2 / g
       endif
    enddo

  endsubroutine calc_kappa_H_ff
!************************************************************************************
  subroutine bremsstrahlung_absorptionCoeff(p,i,frequency,zprotons,kappaRho)

    real :: gaunt,zprotons,ln_physQuant,kappaRho,frequency
    integer :: i
    type (pillar_case) :: p

    intent(in) :: zprotons
    intent(out) :: kappaRho
    
    gaunt = log(exp(5.960 - sqrt3*pi1 * log(frequency*1e-9 * (p%T(i)*1e-4)**(-1.5))) + exp1)

    ln_physQuant = log(p%ne(i)*(p%nHI(i)+p%nHII(i))) + &  !nelectrons*nprotons
         2*log(zprotons) - &
         0.5*log(p%T(i)) - &
         3*log(frequency) + &
         log((1 - exp(-h_planck*frequency*k1_cgs*p%T1(i))))

    kappaRho = bremsstrahlung_constant * exp(ln_physQuant) * gaunt

  endsubroutine bremsstrahlung_absorptionCoeff
!************************************************************************************
  function get_kappa_Hm_bf(p,iw) result(kappa_Hm_bf)
!    """H- ion bound-free continuum opacity.
!    Calculates the sum of the fit function
!    
!        alpha = 1e-18 * sum_m=0^6(a_m * wave^m)
!    
!    Parameters: waves : 1D array
!                    array of wavelengths. Single values
!                    must be contained in an array object.
!    
!    """
    real, dimension(nz) :: kappa_Hm_bf
    integer :: iw
    type (pillar_case) :: p

    ! fit coefficients
    if (wa(iw) .lt. 16000) then
       kappa_Hm_bf = p%hm_bf_factor * p%stim_factor * p%ionization_factor * s_coeff(iw)
    else
       kappa_Hm_bf=0
    endif
!
  endfunction  get_kappa_Hm_bf
!************************************************************************************
  function get_kappa_Hm_ff(p,iw) result(kappa_Hm_ff)
!
!    """H- ion free-free interaction cross section per HI per unit                                 
!    electron pressure.                                                                            
!    """
    real, dimension(nz) :: kappa_Hm_ff,flam
    integer :: iw
    type (pillar_case) :: p
!    
    ! fit coefficients
    
    flam = f_coeff(iw,1) + f_coeff(iw,2)*p%lgtheta + f_coeff(iw,3)*p%lgtheta2 - 26

    kappa_Hm_ff = p%ionization_factor * p%ne * k_cgs * p%T * 10**(flam)

  endfunction get_kappa_Hm_ff
!************************************************************************************
  subroutine gaunt(n, wave_cm, gaunt_factor)
!
!    """
!    Calculate Gaunt Factor
!    Parameters: n : int
!                    principle quantum number
!                waves : ndarray
!                    wavelengths in Angstroms
!    """
!    
    integer :: n
    real :: wave_cm
    real :: lam_n,eps2_n,gaunt_factor
!
    intent(in) :: n,wave_cm
    intent(out) :: gaunt_factor
!
    ! reference wavelength lambda_n
    lam_n = n**2 * hcR1
    
    ! gaunt will be non-zero only where wavelengths are shorter than
    !   lam_n (i.e. the ionizing wavelength)
    if (wave_cm .le. lam_n) then
       eps2_n = ((wave_cm/lam_n) - 1)**2

       ! Gaunt factor. Below experession is for lamda <= lambda_n; zero otherwise
       !121./700=0.17285714285714285
       gaunt_factor = 1 - 0.17285714285714285 * (1-eps2_n) * (n*(1+eps2_n))**(-two_thirds)
    else
       gaunt_factor=0.0
    endif
!
  endsubroutine gaunt
!************************************************************************************
  subroutine grey_parameters(p,sigma_grey)

    real, intent(in) :: sigma_grey
    real, dimension(nz) :: alpha_grey
    type (pillar_case) :: p
    
    p%source_function = sigma_sb*p%T**4
    alpha_grey = 3.68e22 * p%T**(-3.5) *  p%rho
    p%opacity = (alpha_grey+sigma_grey)*p%rho
    p%albedo  = sigma_grey / (alpha_grey + sigma_grey)
    p%opacity1 = 1./p%opacity
    p%opacity2 = p%opacity**2

  endsubroutine grey_parameters
!************************************************************************************
  function get_source_function(T1,iw) result(source_function)

    integer :: i,iw
    real, dimension(nz) :: T1,damping_factor
    real, dimension(nz):: source_function

    !damping_factor_const = h_planck*c_light_cgs*k1_cgs
    damping_factor = damping_factor_const*T1*wcm1(iw)
    do i=1,nz
      if (damping_factor(i) > log_overflow_limit) then
        damping_factor(i) = log_overflow_limit
      endif
    enddo
    
    !source_function_const = 2*h_planck*c_light_cgs**2
    source_function = source_function_const*w1cm5(iw) * 1/(exp(damping_factor)-1)
!    
  endfunction get_source_function
!************************************************************************************
endmodule ContinuousOpacity
