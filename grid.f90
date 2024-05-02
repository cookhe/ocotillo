module grid

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

! Stefan-Boltzmann constant
  real :: sigma_sb=5.670374419d-5
!  
end module grid
