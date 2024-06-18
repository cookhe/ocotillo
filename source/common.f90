module Common

  implicit none

  public

  include 'resolution.inc'
!
! For arrays that need ghost zones
! Number of ghost zones
  integer, parameter :: ng=3

! Number of grid points including ghost zones on each end
  integer, parameter :: mz=nz+2*ng

! Index locations of fist/last zones of computational domain 
  integer, parameter :: n1=ng+1,n2=mz-ng

! Index to store the grey solution
  integer, parameter :: igrey=1
  
! Stefan-Boltzmann constant
  real :: sigma_sb=5.6703744191844314d-05
!  
  real :: float_info_max=3.9085d307

! Mass of proton and inverse
  real :: mp=1.67262192369d-24
  real :: mp1=1 / 1.67262192369d-24

! Mass of electron
  real :: me=9.1093837015d-28

! Planck constant
  real, parameter :: h_planck=6.62607015d-27

! Pi
  real, parameter :: pi=4*atan(1.0)

! Boltzmann constant in erg/K and eV/K
  real :: k_cgs=1.380649d-16,k_eV=8.617333262145179d-05

! Ionization potential energy for Hydrogen in eV
  real :: hydrogen_ionization_eV = 13.6
!
  real, parameter :: c_light_cgs = 2.99792458d10
!
  real :: RydbergEnergy = 2.1798741d-11
!
  real, parameter :: e_electron_cgs  = 4.8032068d-10 ![esu]: charge
!  
  real, parameter :: me_cgs = 9.1093897d-28 ! electron mass cgs
!
  real, parameter :: log10e = log10(exp(1d0))
  real, parameter :: Rcm = 2 * pi**2 * me_cgs * e_electron_cgs**4 / (h_planck**3 * c_light_cgs)
  real, parameter :: Iev = h_planck*c_light_cgs*Rcm * 6.2419d11 ! convert to eV
!  
endmodule Common
