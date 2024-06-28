! Module to calculate the gas state
module GasState

  use Common
  use Columns
  
  implicit none
  private

  public :: calc_hydrogen_ion_frac,solve_gas_state
  public :: get_electron_pressure,read_gas_state_input

  real :: fully_ionized_T
  namelist /gas_state_input/ fully_ionized_T
  
contains
!******************************************
  subroutine read_gas_state_input(inputfile)

    integer :: iu
    character(len=90)   :: inputfile
    
    open(newunit=iu,file=trim(inputfile))
    read(iu,nml=gas_state_input)
    close(iu)

  endsubroutine read_gas_state_input
!******************************************
  subroutine calc_hydrogen_ion_frac(c)
  real, dimension(nz) :: constants,niine_ni,c_ion,exparg
  integer :: i
  type (column_case) :: c
!
! calculate the Saha equation (relative fraction of adjacent ions)
!
  !constants = (sqrt(2*pi*me*k_cgs)*h1_planck)**3   * T**1.5
  !(sqrt(2*pi*me*k_cgs)*h1_planck)**3 = 2414683039571967.0
  constants=2414683039571967. * c%T**1.5
  exparg = c%T1 * hydrogen_ionization_eV*k1_eV
  niine_ni =  constants * exp(-exparg)
    
  c_ion = mp * c%rho1 * niine_ni

  c%NHII_NHINHII = .5*(sqrt(c_ion**2 + 4*c_ion) - c_ion)
!
! Limits due to machine precision mean 4C becomes unresolved next to C**2,
! so manually set the ionization fraction to 1 for temperatures above 20000.
! This is appropriate for our domain where the density is always less than
! ~2e-9 g cm^-3. Would not be appropriate for higher densities.
!
  do i=1,nz
    if (c%T(i) > fully_ionized_T) then 
      c%NHII_NHINHII(i)=1.
    endif
  enddo
!
endsubroutine calc_hydrogen_ion_frac
!******************************************
subroutine solve_gas_state(c)!rho,rho1,NHII_NHINHII,number_density,inv_number_density,&
!     nHI,nHII,ne,ionization_factor)
  !real, dimension(nz), intent(in) :: rho,rho1,NHII_NHINHII
  !real, dimension(nz), intent(out) :: number_density,inv_number_density,nHI,nHII,ne,ionization_factor
  integer :: i
  type (column_case) :: c
  
  c%number_density = c%rho*mp1
  c%inv_number_density = c%rho1*mp
  c%nHII = c%NHII_NHINHII * c%number_density
  c%nHI = c%number_density - c%nHII
  c%ne = c%nHII
  
  do i=1,nz
    if (c%nHI(i) /= 0) then
      c%ionization_factor(i) = c%nHI(i)/(c%nHI(i) + c%nHII(i))
    else ! set to zero.
      c%ionization_factor(i) = 0.
    endif
  enddo

endsubroutine solve_gas_state
!******************************************
function get_electron_pressure(c) result(electron_pressure)
! function  for electron pressure
  !real, intent(in), dimension(nz) :: ne, T
  real, dimension(nz) :: electron_pressure

  type (column_case) :: c
  
  electron_pressure = c%ne * k_cgs * c%T
  
endfunction get_electron_pressure
!******************************************
endmodule GasState
