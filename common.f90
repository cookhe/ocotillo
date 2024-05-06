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
  real :: sigma_sb=5.670374419d-5
!  
  real :: float_info_max=3.9085d307

! Mass of proton and inverse
  real:: mp=1.6726231d-24,mp1=1/mp
!
endmodule common
