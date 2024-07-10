! Module to calculate the gas state
module GasState

  use Common
  use ContinuousOpacity
  
  implicit none
  private

  public :: calc_hydrogen_ion_frac,solve_gas_state
  public :: get_electron_pressure,read_gas_state_input
  public :: wavelength_independent_pillars

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
  subroutine wavelength_independent_pillars(p)
!    
    type (pillar_case) :: p
!    
    p%rho1=1./p%rho
    p%T1=1./p%T
    call calc_hydrogen_ion_frac(p)
    call solve_gas_state(p)
    p%electron_pressure = get_electron_pressure(p)
    p%e_scatter         = get_electron_thomson_scattering(p)
    p%theta             = 5040.* p%T1
    p%theta1            = 1./p%theta
    p%lgtheta           = log10(p%theta)
    p%lgtheta2          = p%lgtheta**2
    p%hm_bf_factor      = get_hydrogen_ion_bound_free(p)
!
  endsubroutine wavelength_independent_pillars
!******************************************
  subroutine calc_hydrogen_ion_frac(p)
  real, dimension(nz) :: constants,niine_ni,c_ion,exparg
  integer :: i
  type (pillar_case) :: p
!
! calculate the Saha equation (relative fraction of adjacent ions)
!
  !constants = (sqrt(2*pi*me*k_cgs)*h1_planck)**3   * T**1.5
  !(sqrt(2*pi*me*k_cgs)*h1_planck)**3 = 2414683039571967.0
  constants=2414683039571967. * p%T**1.5
  exparg = p%T1 * hydrogen_ionization_eV*k1_eV
  niine_ni =  constants * exp(-exparg)
    
  c_ion = mp * p%rho1 * niine_ni

  p%NHII_NHINHII = .5*(sqrt(c_ion**2 + 4*c_ion) - c_ion)
!
! Limits due to machine precision mean 4C becomes unresolved next to C**2,
! so manually set the ionization fraction to 1 for temperatures above 20000.
! This is appropriate for our domain where the density is always less than
! ~2e-9 g cm^-3. Would not be appropriate for higher densities.
!
  do i=1,nz
    if (p%T(i) > fully_ionized_T) then 
      p%NHII_NHINHII(i)=1.
    endif
  enddo
!
endsubroutine calc_hydrogen_ion_frac
!******************************************
subroutine solve_gas_state(p)!rho,rho1,NHII_NHINHII,number_density,inv_number_density,&
!     nHI,nHII,ne,ionization_factor)
  !real, dimension(nz), intent(in) :: rho,rho1,NHII_NHINHII
  !real, dimension(nz), intent(out) :: number_density,inv_number_density,nHI,nHII,ne,ionization_factor
  integer :: i
  type (pillar_case) :: p
  
  p%number_density = p%rho*mp1
  p%inv_number_density = p%rho1*mp
  p%nHII = p%NHII_NHINHII * p%number_density
  p%nHI = p%number_density - p%nHII
  p%ne = p%nHII
  
  do i=1,nz
    if (p%nHI(i) /= 0) then
      p%ionization_factor(i) = p%nHI(i)/(p%nHI(i) + p%nHII(i))
    else ! set to zero.
      p%ionization_factor(i) = 0.
    endif
  enddo

endsubroutine solve_gas_state
!******************************************
function get_electron_pressure(p) result(electron_pressure)
! function  for electron pressure
  !real, intent(in), dimension(nz) :: ne, T
  real, dimension(nz) :: electron_pressure

  type (pillar_case) :: p
  
  electron_pressure = p%ne * k_cgs * p%T
  
endfunction get_electron_pressure
!******************************************
endmodule GasState
