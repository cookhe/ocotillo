! Module to calculate the gas state
module GasState

  use Common
  
  implicit none
  private

  public :: calc_hydrogen_ion_frac,solve_gas_state
  public :: calc_electron_pressure

  real :: fully_ionized_T
  namelist /gas_state_input/ fully_ionized_T
  
contains
!******************************************
subroutine calc_hydrogen_ion_frac(rho,T,NHII_NHINHII)
  real, dimension(nz), intent(in) :: rho,T
  real, dimension(nz), intent(out) :: NHII_NHINHII
  real, dimension(nz) :: constants,niine_ni,C,index
  integer :: i,iu
  
  open(newunit=iu,file='input.in')
  read(iu,nml=gas_state_input)
  close(iu)

  ! calculate the Saha equation (relative fraction of adjacent ions)
  constants = ((sqrt(2*pi*me*k_cgs*T)/h_planck)**3)
  index = hydrogen_ionization_eV/(k_eV*T)
  ! index = 100.
  ! print*, 'E/k_eV*T = ', index
  niine_ni =  constants * exp(-index)
    
  C = mp / rho * niine_ni

  NHII_NHINHII = (sqrt(C**2 + 4*C) - C)/2.

  ! Limits due to machine precision mean 4C becomes unresolved next to C**2,
  ! so manually set the ionization fraction to 1 for temperatures above 20000.
  ! This is appropriate for our domain where the density is always less than
  ! ~2e-9 g cm^-3. Would not be appropriate for higher densities.
  do i=1,nz
    ! print*, T(i), fully_ionized_T
    if (T(i) > fully_ionized_T) then 
      NHII_NHINHII(i)=1.
    endif
  enddo
  ! print*, NHII_NHINHII
  ! NHII_NHINHII = merge(1.0, NHII_NHINHII, T>fully_ionized_T)
  
endsubroutine calc_hydrogen_ion_frac
!******************************************
subroutine solve_gas_state(rho,NHII_NHINHII,n,nHI,nHII,ne,ionization_factor)
  real, dimension(nz), intent(in) :: rho,NHII_NHINHII
  real, dimension(nz), intent(out) :: n,nHI,nHII,ne,ionization_factor
  integer :: i

  n = rho*mp1
  nHII = NHII_NHINHII * n
  nHI = n - nHII
  ne = nHII
  
  do i=1,nz
    if (nHI(i) /= 0) then
      ionization_factor(i) = 1./(1. + nHII(i)/nHI(i))
    else ! set to zero.
      ionization_factor(i) = 0.
    endif
  enddo

endsubroutine solve_gas_state

! Subroune for electron pressure
subroutine calc_electron_pressure(ne, T, electron_pressure)
  real, intent(in), dimension(nz) :: ne, T
  real, intent(out), dimension(nz) :: electron_pressure

  electron_pressure = ne * k_cgs * T
  
endsubroutine calc_electron_pressure

endmodule GasState
