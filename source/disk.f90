! Module for disk model generation/manipulation.
module Disk

  use Common

  implicit none
  private

  public :: read_temperature_input,read_density_input
  public :: calc_grid,calc_density,calc_temperature,calc_wavelength
  
  real :: isoTemp,sigma,midTemp,floorTemp,switchTemp
  real :: rho0,rho_floor,H
  
  namelist /temperature_input/ isoTemp,sigma,midTemp,floorTemp
  namelist /density_input/ rho0,rho_floor,H

contains
!************************************************************************************
   subroutine read_temperature_input(inputfile)
!
    character(len=90)   :: inputfile
!
    open(30,file=trim(inputfile))
    read(30,nml=temperature_input)
    close(30)
!
    switchTemp = (1+epsi) * isoTemp ! must be larger than mdiTemp
!
   endsubroutine read_temperature_input
!************************************************************************************
   subroutine read_density_input(inputfile)
!
    character(len=90)   :: inputfile
!     
    open(30,file=trim(inputfile))
    read(30,nml=density_input)
    close(30)
!
   endsubroutine read_density_input
!************************************************************************************
   subroutine calc_grid(z1,z0,z,dz)
!
      real, dimension(mz) :: z
      real :: dz,z1,z0
      integer :: i
!
      dz = (z1-z0)/(nz-1)
      do i=1,nz
        z(n1+i-1) = z0 + (i-1)*dz
      enddo
!
! Fill in ghost zones
!
      print*,z(n1),z(n2)
!
      do i=1,ng
         z(n1-i) = z(n1) - i*dz
         z(n2+i) = z(n2) + i*dz
      enddo
!
    endsubroutine calc_grid
!************************************************************************************
    subroutine calc_wavelength(w1,w0,wa,wa1,wc,wc1)
!
      real, dimension(nw), intent(out) :: wa,wa1,wc,wc1
      real, intent(in) :: w1,w0
      real :: dw
      integer :: i
!
      if (nw >1 ) then
        dw = (w1-w0)/(nw-1)
        do i=1,nw
          !wavelengths in angstrom
           wa(i) = w0 + (i-1)*dw
        enddo
      else
        wa=w0
      endif
      wa1 = 1./wa
      wc  = wa*1d-8 !wavelengths in cm
      wc1 = 1./wc
      print*, 'waves_angstrom min/max', minval(wa),maxval(wa)
!
    endsubroutine calc_wavelength
!************************************************************************************
   subroutine calc_temperature(T,z,lfrom_read_athena)
!    
    real, dimension(nz), intent(inout) :: T
    real, dimension(mz), intent(in) :: z
    real, dimension(nz), save :: Tgauss
    integer :: i
    logical, optional :: lfrom_read_athena
    logical :: lathena
    logical, save :: lfirst_call=.true.
!
    if (present(lfrom_read_athena)) then
      lathena=lfrom_read_athena
    else
      lathena=.false.
    endif
!
! Check if it's calling from read_athena. Otherwise
! assign to the midTemp value from input.in
!
    if (.not.lathena) T = midTemp
!    
! Calculate the gaussian profile
!
    if (lfirst_call) then
      Tgauss = midTemp * exp(-0.5 * z(n1:n2)**2/sigma**2)
      lfirst_call=.false.
    endif
!        
! If the temperature is greater than the switch temperature,
! it will use the current value (e.g. the Athena input.
! Else choose between the gaussian profile or the floor.
!
    do i=1,nz
      if (T(i) < switchTemp) then
        if (Tgauss(i) > floorTemp) then
          T(i) = Tgauss(i)
        else
          T(i) = floorTemp
        endif
      endif
    enddo
!
    if (.not.lathena) & ! if read from read_athena, it already prints there
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
    rho = rho0*exp(-.5*z(n1:n2)**2/H**2)
!
! Implement a density floor
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
endmodule Disk
