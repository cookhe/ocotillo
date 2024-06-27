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

  real :: switch_ionfraction=1e-2
  real :: ln_const = 19.72694361375364 !log(4./3) + 6*log(e) - log(me*h*c) + 0.5*log(2*pi) - 0.5*log(3*me*kb) 
  real, dimension(7) :: a_coeff
  real, dimension(3,5) :: b_coeff
  real :: AHbf = 1.0449e-26
  real, dimension(6) :: n1_array
  real :: Rangstrom=1.0968e-3 ! Rcm = 2 * pi**2 * me * e**4 / (h**3 * c) 
  real :: alpha0=1.0443e-26     ! absorption per electron
  real, dimension(nw,3) :: f
  
contains
!************************************************************************************
  subroutine pre_calc_opacity_quantities(waves)

    real, dimension(nw) :: waves,lgwave
    integer :: n,iw,i,j

    a_coeff = (/+1.99654,&
         -1.18267e-6,&
         +2.62423e-7,&
         -4.40524e-11,&
         +3.23992e-15,&
         -1.39568e-19,&
         +2.78701e-24/)

    b_coeff = transpose(reshape(                                   &
         (/-2.276300,-1.685000,+0.766610,-0.053356,+0.000000,&
         +15.28270,-9.284600,+1.993810,-0.142631,+0.000000,  &
         -197.789,+190.266,-67.9775,+10.6913,-0.62515/),     &
         (/ size(b_coeff, 2), size(b_coeff, 1) /)))

    do n=1,6
       n1_array(n)=1./n
    enddo
!
    lgwave=log10(waves)
    do iw=1,nw
      do i=1,3
        f(iw,i)=0.
        do j=1,5
          f(iw,i) = f(iw,i) + lgwave(iw)**(j-1) * b_coeff(i,j)
        enddo
      enddo
    enddo
!
  endsubroutine pre_calc_opacity_quantities
!************************************************************************************
  function get_electron_thomson_scattering(inv_number_density, nHII) result(e_scatter)
    real, intent(in), dimension(nz) :: inv_number_density, nHII
    real, dimension(nz) ::  e_scatter
    real :: alpha_e = 0.6648e-24 ! coefficient

    e_scatter = alpha_e * nHII * inv_number_density

  endfunction get_electron_thomson_scattering
!************************************************************************************
  function get_hydrogen_ion_bound_free(electron_pressure,theta) result(hm_bf_factor)

    real, intent(in), dimension(nz) :: electron_pressure,theta
    real, dimension(nz) :: hm_bf_factor
    
    hm_bf_factor=4.158e-10 * electron_pressure * theta**2.5 * 10**(0.754 * theta)

  endfunction get_hydrogen_ion_bound_free
!************************************************************************************
  function get_hydrogen_stimulated_emission(wave1, theta) result(stim_factor)
    real, intent(in) :: wave1  ! MUST be in angstroms
    real, intent(in), dimension(nz) :: theta  ! theta = 5040./T[K]
    real, dimension(nz) :: stim_factor
    real :: chi

    chi = 1.2398e4 * wave1
    stim_factor = 1 - 10**(-chi*theta)

  endfunction get_hydrogen_stimulated_emission
!************************************************************************************
  subroutine calc_opacity_and_albedo(e_scatter,rho,rho1,ne,NHII_NHINHII,nHI,nHII,&
       temp,temp1,theta,theta1,lgtheta,lgtheta2,wave_angstrom,wave_cm,nu_Hz,hm_bf_factor,&
       stim_factor,ionization_factor,opacity,albedo,iw)

    real, dimension(nz) :: e_scatter,rho, rho1,temp, temp1, theta, theta1, lgtheta, lgtheta2
    real, dimension(nz) :: ne, NHII_NHINHII,nHI,nHII
    real, dimension(nz) :: hm_bf_factor, stim_factor, ionization_factor
    real, dimension(nz) :: opacity, albedo
    real :: wave_angstrom, wave_cm, nu_Hz
    integer :: iw

    real, dimension(nz) :: kappa_rad,kappa_H_bf,kappa_Hm_bf,kappa_Hm_ff
    
    intent(in)   :: e_scatter,rho,rho1,ne,temp,temp1,theta,theta1
    intent(in)   :: wave_angstrom,wave_cm,nu_Hz,hm_bf_factor,stim_factor,ionization_factor
    intent(out)  :: opacity, albedo

    kappa_H_bf  = get_kappa_H_bf(wave_angstrom, wave_cm, temp, temp1, stim_factor * ionization_factor)
    kappa_Hm_bf = get_kappa_Hm_bf(wave_angstrom, hm_bf_factor * stim_factor * ionization_factor)
    kappa_Hm_ff = get_kappa_Hm_ff(ne,temp,lgtheta,lgtheta2,ionization_factor,iw)

    kappa_rad = (kappa_H_bf + kappa_Hm_bf + kappa_Hm_ff + e_scatter)*mp1
    opacity = kappa_rad * rho

    call calc_kappa_H_ff(wave_angstrom, nu_Hz, temp, temp1, theta, theta1, stim_factor * ionization_factor,&
         rho,rho1,NHII_NHINHII,ne,nHI,nHII,kappa_rad,opacity)

    albedo = e_scatter * mp1 / kappa_rad

  endsubroutine calc_opacity_and_albedo
!************************************************************************************
  function get_kappa_H_bf(waves_ang, waves_cm, temp, temp1, factor) result(kappa_H_bf)
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
    real, dimension(nz) :: sm, temp, temp1, ktemp, ktemp1, factor, C, kappa_H_bf
    integer :: m=6, n
    real :: chi_n,chi_m,waves_ang,waves_cm,g
!
    intent(in) :: waves_ang, waves_cm, temp, temp1, factor
!
    ! sum for the first m-1 excitation states
    sm = 0.
    ktemp = k_cgs*temp
    ktemp1 = k1_cgs*temp1
    
    do n=1,m
      chi_n = RydbergEnergy * (1. - 1.*n1_array(n)**2)
      call gaunt(n, waves_cm, g)
      sm = sm + g*n1_array(n)**3 * exp(-chi_n*ktemp1)
    enddo
    
    ! excitation energy of state level m
    chi_m = RydbergEnergy * (1. - 1./m**2)
    
    ! Unsold approximation integral
    C = .5*ktemp*RydbergEnergy1 * ( exp(-chi_m*ktemp1) - exp(-RydbergEnergy*ktemp1) )

    kappa_H_bf = AHbf * factor * waves_ang**3 * (C + sm)
    
  endfunction get_kappa_H_bf
!************************************************************************************
  subroutine calc_kappa_H_ff(waves, nu, temp, temp1, theta, theta1, factor,&
       rho,rho1,NHII_NHINHII,ne,nHI,nHII,kappa_rad,opacity)
!
!    """Hydrogen free-free absorption coefficient.
!    
!    Units cm^2 per neutral hydrogen atom.
!    """
    real, dimension(nz) :: kappa_rad,opacity
    real, dimension(nz) :: temp,temp1,theta,theta1,factor,rho,rho1,NHII_NHINHII,ne,nHI,nHII
    real :: g_ff,energyfactor,kappa_H_ff,opacity_bremsstrahlung
    real :: chi1,waves,nu
    integer :: i
!
    intent(in) :: waves,nu,temp, temp1, theta, theta1, factor, rho,rho1,NHII_NHINHII,ne,nHI,nHII
    intent(inout) :: kappa_rad,opacity
!    
    !chi1 = waves/1.2398e4
    chi1 = waves*8.0658170672689143d-5

    do i=1,nz
       if ((1-NHII_NHINHII(i)) .gt. switch_ionfraction) then
          ! Use Gray 2022 function that depends on hydrogen's ionization state.
          g_ff = 1 + 0.3456 * (waves * Rangstrom)**(-one_third) * (theta1(i)*log10e*chi1 + 0.5)
          energyfactor = .5*log10e*theta1(i)*Iev1 * 10**(-theta(i)*Iev)
          kappa_H_ff = factor(i) * alpha0 * waves**3 * g_ff * energyfactor
          kappa_rad(i) = kappa_rad(i) + kappa_H_ff * mp1 ! cm^2 / g
          opacity(i) = opacity(i) + kappa_rad(i) * rho(i) ! 1/cm
       else
          ! Use espression for fully ionized gas.
          !waves_cm = waves / 1e8 ! convert angstroms to centimeters
          !nu_Hz = c*1e8 / waves
          call bremsstrahlung_absorptionCoeff(nu,temp(i),temp1(i),&
               ne(i),nHI(i)+nHII(i),1.0,opacity_bremsstrahlung)
          opacity(i) = opacity(i) + opacity_bremsstrahlung
          kappa_rad(i) = kappa_rad(i) + opacity(i) * rho1(i) ! cm^2 / g
       endif
    enddo

  endsubroutine calc_kappa_H_ff
!************************************************************************************
  subroutine bremsstrahlung_absorptionCoeff(frequency, temperature, temperature1, &
       nelectrons, nprotons, zprotons,kappaRho)

    real :: gaunt
    real :: frequency, temperature, temperature1, nelectrons, nprotons, zprotons
    real :: ln_physQuant,kappaRho

    intent(in) :: frequency, temperature, temperature1, nelectrons, nprotons, zprotons
    intent(out) :: kappaRho
    
    gaunt = log(exp(5.960 - sqrt3*pi1 * log(frequency*1e-9 * (temperature*1e-4)**(-1.5))) + exp1)

    ln_physQuant = log(nelectrons*nprotons) + &
         2*log(zprotons) - &
         0.5*log(temperature) - &
         3*log(frequency) + &
         log((1 - exp(-h_planck*frequency*k1_cgs*temperature1)))

    kappaRho = exp(ln_const) * exp(ln_physQuant) * gaunt

  endsubroutine bremsstrahlung_absorptionCoeff
!************************************************************************************
  function get_kappa_Hm_bf(waves,factor) result(kappa_Hm_bf)
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
    real, dimension(nz) :: factor,kappa_Hm_bf
    real :: s,waves
    integer :: i

    ! fit coefficients
    if (waves .lt. 16000) then
    
    ! loop through each wavelength and perform the summation
       s=0
       do i=1,7
          s = s + a_coeff(i)*waves**(i-1)
       enddo

       kappa_Hm_bf = factor*1e-17 * s
    else
       kappa_Hm_bf=0
    endif
!
  endfunction  get_kappa_Hm_bf
!************************************************************************************
  function get_kappa_Hm_ff(ne, temp, lgtheta, lgtheta2, factor,iw) result(kappa_Hm_ff)
!
!    """H- ion free-free interaction cross section per HI per unit                                 
!    electron pressure.                                                                            
!    """
    real, dimension(nz) :: ne,temp,kappa_Hm_ff,lgtheta, lgtheta2, factor,flam
    integer :: iw
!    
    intent(in) :: ne,temp, lgtheta, lgtheta2, iw
!
    ! fit coefficients
    
    flam = f(iw,1) + f(iw,2)*lgtheta + f(iw,3)*lgtheta2 - 26

    kappa_Hm_ff = factor * ne * k_cgs * temp * 10**(flam)

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
  subroutine grey_parameters(rho,T,sigma_grey,&
       B_grey,kappa_grey,omega_grey)

    real, intent(in), dimension(nz) :: rho,T
    real, intent(in) :: sigma_grey
    real, dimension(nz) :: alpha_grey
    real, intent(out), dimension(nz) :: B_grey,kappa_grey,omega_grey
    
    B_grey = sigma_sb*T**4
    alpha_grey = 3.68e22 * T**(-3.5) *  rho
    kappa_grey = (alpha_grey+sigma_grey)*rho
    omega_grey = sigma_grey / (alpha_grey + sigma_grey)

  endsubroutine grey_parameters
!************************************************************************************
endmodule ContinuousOpacity
