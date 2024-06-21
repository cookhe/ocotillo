! Module to calculate the gas state
module GasState

  use Common
  
  implicit none
  private

  public :: calc_hydrogen_ion_frac,solve_gas_state
  public :: get_electron_pressure,read_gas_state_input

  real :: fully_ionized_T
  namelist /gas_state_input/ fully_ionized_T
  
contains
!******************************************
  subroutine read_gas_state_input()

    integer :: iu

    open(newunit=iu,file='./input.in')
    read(iu,nml=gas_state_input)
    close(iu)

  endsubroutine read_gas_state_input
!******************************************
  subroutine calc_hydrogen_ion_frac(rho1,T,T1,NHII_NHINHII)
  real, dimension(nz), intent(in) :: rho1,T,T1
  real, dimension(nz), intent(out) :: NHII_NHINHII
  real, dimension(nz) :: constants,niine_ni,C,exparg
  integer :: i
!
! calculate the Saha equation (relative fraction of adjacent ions)
!
  constants = (sqrt(2*pi*me*k_cgs*T)*h1_planck)**3
  exparg = T1 * hydrogen_ionization_eV*k1_eV
  niine_ni =  constants * exp(-exparg)
    
  C = mp * rho1 * niine_ni

  NHII_NHINHII = .5*(sqrt(C**2 + 4*C) - C)
!
! Limits due to machine precision mean 4C becomes unresolved next to C**2,
! so manually set the ionization fraction to 1 for temperatures above 20000.
! This is appropriate for our domain where the density is always less than
! ~2e-9 g cm^-3. Would not be appropriate for higher densities.
!
  do i=1,nz
    if (T(i) > fully_ionized_T) then 
      NHII_NHINHII(i)=1.
    endif
  enddo
!
endsubroutine calc_hydrogen_ion_frac
!******************************************
subroutine solve_gas_state(rho,NHII_NHINHII,number_density,nHI,nHII,ne,ionization_factor)
  real, dimension(nz), intent(in) :: rho,NHII_NHINHII
  real, dimension(nz), intent(out) :: number_density,nHI,nHII,ne,ionization_factor
  integer :: i

  number_density = rho*mp1
  nHII = NHII_NHINHII * number_density
  nHI = number_density - nHII
  ne = nHII
  
  do i=1,nz
    if (nHI(i) /= 0) then
      ionization_factor(i) = 1./(1. + nHII(i)/nHI(i))
    else ! set to zero.
      ionization_factor(i) = 0.
    endif
  enddo

endsubroutine solve_gas_state
!******************************************
function get_electron_pressure(ne, T) result(electron_pressure)
! function  for electron pressure
  real, intent(in), dimension(nz) :: ne, T
  real, dimension(nz) :: electron_pressure

  electron_pressure = ne * k_cgs * T
  
endfunction get_electron_pressure
!******************************************
endmodule GasState
