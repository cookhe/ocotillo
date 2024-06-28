module Common

  implicit none

  public

  include '../resolution.in'

  include './pillars.inc'
!
! For arrays that need ghost zones
! Number of ghost zones
  integer, parameter :: ng=3

! Number of grid points including ghost zones on each end
  integer, parameter :: mz=nz+2*ng

! Index locations of fist/last zones of computational domain 
  integer, parameter :: n1=ng+1,n2=mz-ng

! Stefan-Boltzmann constant
  real :: sigma_sb=5.6703744191844314d-05
!  
  real :: epsi=epsilon(1.0) !machine precision
! huge(1.0) =3.9085d307  
  real :: log_overflow_limit = floor(log10(huge(1.0)))
!
! Mass of proton and inverse
  real, parameter :: mp=1.67262192369d-24
  real :: mp1=1./mp
  real :: amu=1.66053906660e-24

! Mass of electron
  real :: me=9.1093837015d-28

! Planck constant
  real, parameter :: h_planck=6.62607015d-27
  real :: h1_planck = 1./h_planck

! Pi
  real, parameter :: pi=4*atan(1.0)
  real :: pi1=1d0/pi
  real :: sqrt3 = sqrt(3d0)
  real :: exp1 = exp(1d0)

! Boltzmann constant in erg/K and eV/K
  real, parameter :: k_cgs=1.380649d-16
  real, parameter :: k_eV=8.617333262145179d-05
  real :: k1_cgs = 1./k_cgs
  real :: k1_eV = 1./k_eV

! Ionization potential energy for Hydrogen in eV
  real :: hydrogen_ionization_eV = 13.6
!
  real, parameter :: c_light_cgs = 2.99792458d10
!
  real, parameter :: RydbergEnergy = 2.1798741d-11
  real :: RydbergEnergy1 = 1./RydbergEnergy
!
  real, parameter :: e_electron_cgs  = 4.8032068d-10 ![esu]: charge
!  
  real, parameter :: me_cgs = 9.1093897d-28 ! electron mass cgs
!
  real, parameter :: log10e = log10(exp(1d0))
  real, parameter :: Rcm = 2 * pi**2 * me_cgs * e_electron_cgs**4 / (h_planck**3 * c_light_cgs)
  real, parameter :: Iev = h_planck*c_light_cgs*Rcm * 6.2419d11 ! convert to eV
  real :: Iev1 = 1./Iev
  real :: bremsstrahlung_constant = 369234910.67735046 ! exp(log(4./3) + 6*log(e) - log(me*h*c) + 0.5*log(2*pi) - 0.5*log(3*me*kb))
!  
  logical :: lfirst,lroot
!
  real, parameter :: G_Newton_cgs=6.6742d-8
  real :: SolarMass=1.988409870698051e+33
!
  real :: hcR1 =  h_planck*c_light_cgs/RydbergEnergy
  real :: one_third = 1d0/3
  real :: two_thirds = 2d0/3
!
endmodule Common
