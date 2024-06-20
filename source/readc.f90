module ReadAthena

  use Common
  use, intrinsic :: iso_c_binding, only: c_float, c_double  
  implicit none

  private

  public :: read_from_athena

  character(len=90)   :: FileName

  real :: Mbh_SolarMasses,r0ref_rg
  real :: aspect_ratio,mean_molecular_weight,rho0
  
  namelist /athena_input/ FileName,Mbh_SolarMasses,r0ref_rg,&
       aspect_ratio,mean_molecular_weight,rho0

  
contains
!************************************************************************************
    subroutine read_from_athena(z,dz,rhoc,temp)
  
      integer :: coordsys
      integer :: nxc,nyc,nzc

      real, intent(out) :: dz
      real, dimension(mz), intent(out) :: z
      real, dimension(nz,ny,nx), intent(out) :: rhoc,temp

  !order: 
  !coordsys,nx,ny,nz,nvar,nscalars,selfgrav_boolean, particles_boolean,gamma1,cs,t,dt,x,y,z,rho,rux,ruy,ruz,eng
  !
      open(40,file='./input.in')
      read(40,nml=athena_input)
      close(40)

      print*,trim(FileName)
      open(99, file = trim(FileName), form = 'unformatted', ACCESS='stream')
      read(99) coordsys
      read(99) nxc,nyc,nzc
      close(99)
!
      print*,"Athena input x dimensionality=",nxc
      if (nxc /= nx) then
         print*,"RT dimensionality from resolution.in =",nx
         print*,"The Athena x dimensionality does not match the chosen for the RT post-processing. Fix."
         stop
      endif
!
      print*,"Athena input y dimensionality=",nyc
      if (nyc /= ny) then
         print*,"RT dimensionality from resolution.in =",ny
         print*,"The Athena y dimensionality does not match the chosen for the RT post-processing. Fix."
         stop
      endif
!
      print*,"Athena input z dimensionality=",nzc
      if (nzc /= nz) then
         print*,"RT dimensionality from resolution.in =",nz
         print*,"The Athena z dimensionality does not match the chosen for the RT post-processing. Fix."
         stop
      endif
!  
      call read_rest_of_file(z,dz,rhoc,temp,nxc,nyc,nzc,FileName)
!  
    endsubroutine read_from_athena
!************************************************************************************
    subroutine read_rest_of_file(z,dz,rhoc,temp,nxc,nyc,nzc,FileName)
  
      integer :: coordsys
      integer :: nxc,nyc,nzc
      integer :: nvar,nscalars
      integer :: selfgrav_boolean, particles_boolean
      real(c_double)       :: gamma1,csc
      real(c_double)       :: tc,dtc
      real(c_double), dimension(nxc) :: xc
      real(c_double), dimension(nyc) :: yc
      real(c_double), dimension(nzc) :: zc
      real(c_double), dimension(nzc,nyc,nxc) :: ruxc,ruyc,ruzc,engc
      character(len=90)   :: FileName
  
      real, intent(out) :: dz
      real, dimension(mz), intent(out) :: z
      real, dimension(nz,ny,nx), intent(out) :: rhoc,temp
  !order: 
  !coordsys,nx,ny,nz,nvar,nscalars,selfgrav_boolean, particles_boolean,gamma1,cs,t,dt,x,y,z,rho,rux,ruy,ruz,eng
  !
  
      open(99, file = trim(FileName), form = 'unformatted', ACCESS='stream')
      read(99) coordsys
      read(99) nxc,nyc,nzc,nvar,nscalars
      read(99) selfgrav_boolean, particles_boolean
      read(99) gamma1,csc
      read(99) tc,dtc
      read(99) xc
      read(99) yc
      read(99) zc
      read(99) rhoc
      read(99) ruxc
      read(99) ruyc
      read(99) ruzc
      read(99) engc
      close(99)
!
      call postprocess_athena_values(zc,gamma1,rhoc,ruxc,ruyc,ruzc,engc,nxc,nyc,nzc,z,temp,dz)
!
    endsubroutine read_rest_of_file
!************************************************************************************
    subroutine postprocess_athena_values(zc,gamma1,rhoc,ruxc,ruyc,ruzc,engc,nxc,nyc,nzc,z,temp,dz)

      integer :: i,ix,iy
      integer :: nxc,nyc,nzc
      real :: gamma,gamma_inv,rg,Mbh,rr,g0,Omega,H
      real :: unit_time,unit_length,unit_velocity
      real :: unit_density,unit_mass,unit_energy,unit_edens
      real, intent(in) :: gamma1
      real, intent(out) :: dz
      real, dimension(nzc) :: zc
      real, dimension(nz) :: ekin,eint,cs2,rho1
      real, dimension(mz), intent(out) :: z
      real, dimension(nzc,nyc,nxc) :: ruxc,ruyc,ruzc,engc
      real, dimension(nz,ny,nx), intent(inout) :: rhoc
      real, dimension(nz,ny,nx), intent(out) :: temp
!
      dz=(zc(nz)-zc(1))/(nz-1)
      z(n1:n2)=zc
!
! Fill in ghost zones
!
      do i=1,ng
         z(n1-i) = z(n1) - i*dz
         z(n2+i) = z(n2) + i*dz
      enddo
!
      Mbh   = Mbh_SolarMasses*SolarMass
      rg    = 2 * (G_Newton_cgs/c_light_cgs**2) * Mbh
      rr    = r0ref_rg * rg
      g0    = G_Newton_cgs*Mbh
      Omega = sqrt(g0)*rr**(-1.5)
      H = aspect_ratio * rr
!
! code units for conversion
!
      unit_time     = 1./Omega
      unit_length   = H
      unit_velocity = unit_length/unit_time
      unit_density  = rho0
      unit_mass     = unit_density*unit_length**3
      unit_energy   = unit_mass*unit_velocity**2
      unit_edens    = unit_energy/unit_length**3
!
      zc=zc*unit_length
!
      rhoc = rhoc*unit_density
!
      ruxc = ruxc*unit_density*unit_velocity
      ruyc = ruyc*unit_density*unit_velocity
      ruzc = ruzc*unit_density*unit_velocity
!
      gamma = 1+gamma1
      gamma_inv = 1./gamma
!
      engc = engc*unit_energy/unit_length**3
      do ix=1,nx
        do iy=1,ny
          rho1=1./rhoc(:,iy,ix)
          ekin = 0.5*rho1*(ruxc(1:nz,iy,ix)**2 + &
                           ruyc(1:nz,iy,ix)**2 + &
                           ruzc(1:nz,iy,ix)**2)
          eint = engc(:,iy,ix) - ekin
          cs2 = gamma*gamma1 * eint * rho1
          temp(:,iy,ix) = mean_molecular_weight * amu * cs2 * gamma_inv * k1_cgs
        enddo
      enddo
!
      print*,'unit_time     = ',unit_time,' s'
      print*,'unit_length   = ',unit_length,' cm'
      print*,'unit_velocity = ',unit_velocity,' cm/s'
      print*,'unit_density  = ',unit_density,' g/cm3'
      print*,'unit_mass     = ',unit_mass,' g'
      print*,'unit_energy   = ',unit_energy,' erg'
      print*,'unit_edens    = ',unit_edens,' erg/cm3'
!
      print*,'minval(rho),maxval(rho),mean(rho)',&
           minval(rhoc),maxval(rhoc),sum(rhoc)/(nxc*nyc*nzc)
      print*,'minval(T),maxval(T),mean(T)',&
           minval(temp),maxval(temp),sum(temp)/(nxc*nyc*nzc)
!
    endsubroutine postprocess_athena_values
!************************************************************************************
endmodule ReadAthena
