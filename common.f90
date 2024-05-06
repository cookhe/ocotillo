module common

  implicit none

  public

! Number of wavelengths
  integer, parameter :: nw=1

! Number of grid points in computational domain
  integer, parameter :: nz=640
    
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
  real :: sigma_sb=5.6703744191844314e-05
!  
  real :: float_info_max=3.9085d307

! Mass of proton and inverse
  real :: mp=1.67262192369e-24
  real :: mp1=1 / 1.67262192369e-24

! Mass of electron
  real :: me=9.1093837015e-28

! Planck constant
  real :: h_planck=6.62607015e-27

! Pi
  real :: pi=4*atan(1.0)

! Boltzmann constant in erg/K and eV/K
  real :: k_cgs=1.380649e-16,k_eV=8.617333262145179e-05

! Ionization potential energy for Hydrogen in eV
  real :: hydrogen_ionization_eV = 13.6
!
endmodule common
