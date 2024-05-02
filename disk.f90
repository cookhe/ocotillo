! Module for disk model generation/manipulation.
module disk

  use grid

  implicit none
  private

  public :: calc_grid,calc_density,calc_temperature
  
  real :: eps,isoTemp,sigma,midTemp,floorTemp,switchTemp
  real :: rho0,rho_floor,H
  
  namelist /temperature_input/ eps,isoTemp,sigma,midTemp,floorTemp
  namelist /density_input/ rho0,rho_floor,H
  
contains
!************************************************************************************
    subroutine calc_grid(z1,z0,z,dz)
!
      real, dimension(mz) :: z
      real :: dz,z1,z0
      integer :: i
!
      dz = (z1-z0)/nz
      do i=1,nz
        z(n1+i-1) = z0 + (i-1)*dz
      enddo
!
   endsubroutine calc_grid  
!************************************************************************************
  subroutine calc_temperature(T,z)
!    
    real, dimension(nz), intent(inout) :: T
    real, dimension(mz), intent(in) :: z
    real, dimension(nz) :: Tgauss
    integer :: i
!
    open(30,file='input.in')
    read(30,nml=temperature_input)
    close(30)
!
    switchTemp = (1+eps) * isoTemp ! must be larger than mdiTemp
!
! In the future, this T will be the Athena output
!
    T = midTemp    
!    
! Calculate the gaussian profile
!
    Tgauss = midTemp * exp(-0.5 * z(n1:n2)**2/sigma**2)
!        
! Choose between the gaussian profile, the floor, or the current value.
!
    do i=1,nz
      if ((T(i) < switchTemp) .and. (Tgauss(i) > floorTemp)) then
        T(i) = Tgauss(i)
      elseif ((T(i) < switchTemp) .and. (Tgauss(i) < floorTemp)) then
        T(i) = floorTemp
      endif
    enddo
!
    print*,"maxval(T), minval(T)",maxval(T), minval(T)
!
    endsubroutine calc_temperature
!************************************************************************************
   subroutine calc_density(rho,z)
!
    real, dimension(nz), intent(inout) :: rho
    real, dimension(mz), intent(in) :: z
    integer :: i
!
    open(40,file='input.in')
    read(40,nml=density_input)
    close(40)
!
    rho = rho0*exp(-.5*z(n1:n2)**2/H**2)
!
! Implement a desity floor    
!
    do i=1,nz
      if (rho(i) < rho_floor) then 
        rho(i) = rho_floor
      endif
    enddo
!
    print*,"maxval(rho), minval(rho)",maxval(rho), minval(rho)
!    
  endsubroutine calc_density
!************************************************************************************
endmodule disk
