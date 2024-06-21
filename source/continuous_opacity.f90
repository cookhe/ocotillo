! Module to calculate hydrogen cross sections as described
!   in Gray's (2022) description of continuous opacity in
!   stellar atmospheres.
module ContinuousOpacity

  use Common

  implicit none
  private

  public :: get_hydrogen_stimulated_emission
  public :: get_electron_thomson_scattering,get_theta
  public :: get_hydrogen_ion_bound_free
  public :: calc_opacity_and_albedo
  public :: grey_parameters

contains
!************************************************************************************
  function get_electron_thomson_scattering(number_density, nHII) result(e_scatter)
    real, intent(in), dimension(nz) :: number_density, nHII
    real, dimension(nz) ::  e_scatter
    real :: alpha_e = 0.6648e-24 ! coefficient

    e_scatter = alpha_e * (nHII / number_density)  

  endfunction get_electron_thomson_scattering
!************************************************************************************
  function get_theta(T1) result(theta)
    real, intent(in), dimension(nz) :: T1
    real, dimension(nz) ::  theta
    
    theta = 5040.* T1
    
  endfunction get_theta
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
       temp,temp1,theta,wave_angstrom,hm_bf_factor,stim_factor,ionization_factor,opacity,albedo)

    real, dimension(nz) :: e_scatter,rho, rho1,temp, temp1, theta, ne, NHII_NHINHII,nHI,nHII
    real, dimension(nz) :: hm_bf_factor, stim_factor, ionization_factor
    real, dimension(nz) :: opacity, albedo
    real :: wave_angstrom

    real, dimension(nz) :: kappa_rad,kappa_H_bf,kappa_Hm_bf,kappa_Hm_ff
    
    intent(in)   :: e_scatter,rho,rho1,ne,temp,temp1,wave_angstrom,hm_bf_factor,stim_factor,ionization_factor
    intent(out)  :: opacity, albedo

    kappa_H_bf  = get_kappa_H_bf(wave_angstrom, temp, temp1, stim_factor * ionization_factor)
    kappa_Hm_bf = get_kappa_Hm_bf(wave_angstrom, hm_bf_factor * stim_factor * ionization_factor)
    kappa_Hm_ff = get_kappa_Hm_ff(ne,wave_angstrom,temp,theta,ionization_factor)

    kappa_rad = (kappa_H_bf + kappa_Hm_bf + kappa_Hm_ff + e_scatter)*mp1
    opacity = kappa_rad * rho

    call calc_kappa_H_ff(wave_angstrom, temp, theta, stim_factor * ionization_factor,&
         rho,rho1,NHII_NHINHII,ne,nHI,nHII,kappa_rad,opacity)

    albedo = e_scatter / (kappa_rad*mp)

  endsubroutine calc_opacity_and_albedo
!************************************************************************************
  function get_kappa_H_bf(waves, temp, temp1, factor) result(kappa_H_bf)
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
    real :: A,chi_n,chi_m,waves, g, n1
!
    intent(in) :: waves, temp, temp1, factor
!
    A = 1.0449e-26    ! cm^2 A^-3 - fundamental constants
           
    ! sum for the first m-1 excitation states
    sm = 0.
    ktemp = k_cgs*temp
    ktemp1 = k1_cgs*temp1
    
    do n=1,m
      n1=1./n
      chi_n = RydbergEnergy * (1. - 1.*n1**2)
      call gaunt(n, waves, g)
      sm = sm + g*n1**3 * exp(-chi_n*ktemp1)
    enddo
    
    ! excitation energy of state level m
    chi_m = RydbergEnergy * (1. - 1./m**2)
    
    ! Unsold approximation integral
    C = .5*ktemp*RydbergEnergy1 * ( exp(-chi_m*ktemp1) - exp(-RydbergEnergy*ktemp1) )

    kappa_H_bf = A * factor * waves**3 * (C + sm)
    
  endfunction get_kappa_H_bf
!************************************************************************************
  subroutine calc_kappa_H_ff(waves, temp, theta, factor,&
    rho,rho1,NHII_NHINHII,ne,nHI,nHII,&
    kappa_rad,opacity)
!
!    """Hydrogen free-free absorption coefficient.
!    
!    Units cm^2 per neutral hydrogen atom.
!    """
    real, dimension(nz) :: kappa_rad,opacity
    real, dimension(nz) :: temp,theta,theta1,factor,rho,rho1,NHII_NHINHII,ne,nHI,nHII
    real :: g_ff,energyfactor,kappa_H_ff,opacity_bremsstrahlung
    real :: alpha0=1.0443e-26     ! absorption per electron
    real :: Rangstrom=1.0968e-3   ! Rcm = 2 * pi**2 * me * e**4 / (h**3 * c)
    real :: chi1,waves,switch_ionfraction
    integer :: i
!
    intent(in) :: waves, temp, factor, rho,rho1,NHII_NHINHII,ne,nHI,nHII
    intent(inout) :: kappa_rad,opacity
!    
    switch_ionfraction=1e-2

    chi1 = waves/1.2398e4
    theta1 = temp/5040

    do i=1,nz
       if ((1-NHII_NHINHII(i)) .gt. switch_ionfraction) then
          ! Use Gray 2022 function that depends on hydrogen's ionization state.
          g_ff = 1 + 0.3456 * (waves * Rangstrom)**(-1./3) * (theta1(i)*log10e*chi1 + 0.5)
          energyfactor = .5*log10e*theta1(i)*Iev1 * 10**(-theta(i)*Iev)
          kappa_H_ff = factor(i) * alpha0 * waves**3 * g_ff * energyfactor
          kappa_rad(i) = kappa_rad(i) + kappa_H_ff * mp1 ! cm^2 / g
          opacity(i) = opacity(i) + kappa_rad(i) * rho(i) ! 1/cm
       else
          ! Use espression for fully ionized gas.
          !waves_cm = waves / 1e8 ! convert angstroms to centimeters
          !nu_Hz = c*1e8 / waves
          call bremsstrahlung_absorptionCoeff(c_light_cgs*1e8/waves,temp(i),ne(i),nHI(i)+nHII(i),1.0,opacity_bremsstrahlung)
          opacity(i) = opacity(i) + opacity_bremsstrahlung
          kappa_rad(i) = kappa_rad(i) + opacity(i) * rho1(i) ! cm^2 / g
       endif
    enddo

  endsubroutine calc_kappa_H_ff
!************************************************************************************
  subroutine bremsstrahlung_absorptionCoeff(frequency, temperature, nelectrons, nprotons, zprotons,kappaRho)

    real :: h,kb,ln_const,gaunt
    real :: frequency, temperature, nelectrons, nprotons, zprotons
    real :: ln_physQuant,kappaRho

    intent(in) :: frequency, temperature, nelectrons, nprotons, zprotons
    intent(out) :: kappaRho
    
    h=h_planck
    kb=k_cgs
    
    ln_const = 19.72694361375364 !log(4./3) + 6*log(e) - log(me*h*c) + 0.5*log(2*pi) - 0.5*log(3*me*kb)

    gaunt = log(exp(5.960 - sqrt(3.)/pi * log(frequency/1e9 * (temperature/1e4)**(-1.5))) + exp(1.))

    ln_physQuant = log(nelectrons*nprotons) + &
         2*log(zprotons) - &
         0.5*log(temperature) - &
         3*log(frequency) + &
         log((1 - exp(-h*frequency/(kb*temperature))))

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
    real, dimension(7) :: a
    real :: s,waves
    integer :: i

    ! fit coefficients
    if (waves .lt. 16000) then
       a = (/+1.99654,&
       -1.18267e-6,&
       +2.62423e-7,&
       -4.40524e-11,&
       +3.23992e-15,&
       -1.39568e-19,&
       +2.78701e-24/)
    
    ! loop through each wavelength and perform the summation
       s=0
       do i=1,7
          s = s + a(i)*waves**(i-1)
       enddo

       kappa_Hm_bf = factor*1e-17 * s
    else
       kappa_Hm_bf=0
    endif
!
  endfunction  get_kappa_Hm_bf
!************************************************************************************
  function get_kappa_Hm_ff(ne, waves, temp, theta, factor) result(kappa_Hm_ff)
!
!    """H- ion free-free interaction cross section per HI per unit                                 
!    electron pressure.                                                                            
!    """
    real, dimension(nz) :: ne,temp,kappa_Hm_ff,theta,factor,flam
    real, dimension(3,5) :: b
    real, dimension(3) :: f
    real :: waves
    integer i,j
!    
    intent(in) :: ne,waves,temp, theta
!
    ! fit coefficients
    
    b = transpose(reshape(                                   &
         (/-2.276300,-1.685000,+0.766610,-0.053356,+0.000000,&
         +15.28270,-9.284600,+1.993810,-0.142631,+0.000000,  &
         -197.789,+190.266,-67.9775,+10.6913,-0.62515/),     &
         (/ size(b, 2), size(b, 1) /)))

    do i=1,3
      f(i)=0.
      do j=1,5
        f(i) = f(i) + log10(waves)**(j-1) * b(i,j)
      enddo
    enddo

    flam = f(1) + f(2)*log10(theta) + f(3)*log10(theta)**2 - 26

    kappa_Hm_ff = factor * ne * k_cgs * temp * 10**(flam)

  endfunction get_kappa_Hm_ff
!************************************************************************************
  subroutine gaunt(n, waves, gaunt_factor)
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
    real :: waves, waves_cm
    real :: R,c,h,lam_n,eps_n,gaunt_factor
!
    intent(in) :: n,waves
    intent(out) :: gaunt_factor
!
    R = RydbergEnergy
    c = c_light_cgs
    h = h_planck
    
    ! convert Angstroms to centimeters
    waves_cm = waves*1e-8
    
    ! reference wavelength lambda_n
    lam_n = n**2 * h*c/R
    
    ! gaunt will be non-zero only where wavelengths are shorter than
    !   lam_n (i.e. the ionizing wavelength)
    if (waves_cm .le. lam_n) then
       eps_n = (waves_cm/lam_n) - 1

       ! Gaunt factor. Below experession is for lamda <= lambda_n; zero otherwise
       gaunt_factor = 1 - (121./700.) * (1-eps_n**2) / (n*(1+eps_n**2))**(2./3)       
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
