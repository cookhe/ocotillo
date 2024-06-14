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
  subroutine calc_opacity_and_albedo(e_scatter,rho,ne,temp,wave_angstrom,&
       hm_bf_factor,stim_factor,ionization_factor,opacity,albedo)

    real, dimension(nz) :: e_scatter,rho, temp, ne
    real, dimension(nz) :: hm_bf_factor, stim_factor, ionization_factor
    real, dimension(nz) :: opacity, albedo
    real :: wave_angstrom

    real, dimension(nz) :: kappa_rad,kappa_H_bf,kappa_H_ff,kappa_Hm_bf,kappa_Hm_ff
    
    intent(in)   :: e_scatter,rho,ne,temp,wave_angstrom,hm_bf_factor,stim_factor,ionization_factor
    intent(out)  :: opacity, albedo

    call calc_kappa_H_bf (wave_angstrom, temp, stim_factor * ionization_factor,  kappa_H_bf)
    call calc_kappa_H_ff (wave_angstrom, temp, stim_factor * ionization_factor,  kappa_H_ff)
    call calc_kappa_Hm_bf(wave_angstrom, hm_bf_factor * stim_factor * ionization_factor, kappa_Hm_bf)
    call calc_kappa_Hm_ff(ne,wave_angstrom,temp,ionization_factor,kappa_Hm_ff)
    
    kappa_rad = kappa_H_bf + kappa_H_ff + kappa_Hm_bf + kappa_Hm_ff + e_scatter

    albedo = e_scatter / kappa_rad
    opacity = kappa_rad * rho / mp

  endsubroutine calc_opacity_and_albedo
!******************************************
  subroutine calc_kappa_H_bf(waves, temp, factor, kappa_H_bf)
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
    real, dimension(nz) :: sm, temp, ktemp, ktemp1, factor, C, kappa_H_bf
    integer :: m=6, n
    real :: A,R,chi_n,chi_m,waves, g
!
    intent(in) :: waves, temp, factor
    intent(out) :: kappa_H_bf
!
    A = 1.0449e-26    ! cm^2 A^-3 - fundamental constants
    R = RydbergEnergy ! erg       - Rydberg energy 2.1798741e-11
           
    ! sum for the first m-1 excitation states
    sm = 0.
    ktemp = k_cgs*temp
    ktemp1 = 1./ktemp
    
    do n=1,m
      chi_n = R * (1. - 1./n**2)
      call gaunt(n, waves, g)
      sm = sm + g/n**3 * exp(-chi_n*ktemp1)
     enddo
    
    ! excitation energy of state level m
    chi_m = R * (1. - 1./m**2)
    
    ! Unsold approximation integral
    C = ktemp/(2*R) * ( exp(-chi_m*ktemp1) - exp(-R*ktemp1) )
    
    kappa_H_bf = A * factor * waves**3 * (C + sm)
    
  endsubroutine calc_kappa_H_bf
!******************************************
  subroutine calc_kappa_H_ff(waves, temp, factor, kappa_H_ff)
!    """Hydrogen free-free absorption coefficient.
!    
!    Units cm^2 per neutral hydrogen atom.
!    """
    real, dimension(nz) :: temp,theta,g_ff,energyfactor,factor,kappa_H_ff
    real :: I,k,R,c,h,e,me,alpha0,Rangstrom,chi,log10e,waves
!
    intent(in) :: waves, temp, factor
    intent(out) :: kappa_H_ff
!    
    k = k_cgs
    R = RydbergEnergy ! erg       - Rydberg energy 2.1798741e-11
    c = c_light_cgs
    h = h_planck
    e  = e_electron_cgs
    me = me_cgs
    I = Iev
           
    ! absorption per electron
    alpha0 = 1.0443e-26
    !Rcm = 2 * pi**2 * me * e**4 / (h**3 * c)
    Rangstrom = 1.0968e-3
    ! gaunt factor
    theta = 5040/temp
    chi = 1.2398e4 / waves

    g_ff = 1 + 0.3456 / (waves * Rangstrom)**(1./3) * (log10e/(theta*chi) + 0.5)
    
    energyfactor = log10e/(2*theta*I) * 10**(-theta*I)
    
    kappa_H_ff = factor * alpha0 * waves**3 * g_ff * energyfactor

  endsubroutine calc_kappa_H_ff
!******************************************
  subroutine calc_kappa_Hm_bf(waves,factor,kappa_Hm_bf)
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
    if (waves > 16000) then 
       a = (/+1.99654,&
       -1.18267e-6,&
       +2.62423e-7,&
       -4.40524e-11,&
       +3.23992e-15,&
       -1.39568e-19,&
       +2.78701e-24/)
    
    ! loop through each wavelength and perform the summation
       s=0
       do i=0,6
          s = s + a(i+1)*waves**i
       enddo

       kappa_Hm_bf = factor*1e-17 * s
    else
       kappa_Hm_bf=0
    endif
!
  endsubroutine  calc_kappa_Hm_bf
!******************************************  
  subroutine calc_kappa_Hm_ff(ne, waves, temp, factor, kappa_Hm_ff)
!
!    """H- ion free-free interaction cross section per HI per unit                                 
!    electron pressure.                                                                            
!    """
    real, dimension(nz) :: ne,temp,kappa_Hm_ff,theta,factor,flam
    real, dimension(3,5) :: b
    real, dimension(3) :: f
    real :: k,waves
    integer i,j
!    
    intent(in) :: ne,waves,temp
    intent(out) :: kappa_Hm_ff
!
    k = k_cgs  ! erg K^-1  - Boltzmann constant                                            
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

    theta = 5040./temp

    flam = f(1) + f(2)*log10(theta) + f(3)*log10(theta)**2 - 26

    kappa_Hm_ff = factor * ne * k * temp * 10**(flam)

  endsubroutine calc_kappa_Hm_ff
!******************************************  
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
!******************************************
endmodule ContinuousOpacity
