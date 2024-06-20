module ReadC

  use Common
  use, intrinsic :: iso_c_binding, only: c_float, c_double  
  implicit none

  private

  public :: read_from_athena

  character(len=90)   :: FileName
  
  namelist /athena_input/ FileName
  
contains
!************************************************************************************
    subroutine read_from_athena(z,dz,rhoc)
  
      integer :: coordsys
      integer :: nxc,nyc,nzc

      real, intent(out) :: dz
      real, dimension(mz), intent(out) :: z
      real, dimension(nz,ny,nx), intent(out) :: rhoc

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
      call read_rest_of_file(z,dz,rhoc,nxc,nyc,nzc,FileName)
!  
    endsubroutine read_from_athena
!************************************************************************************
    subroutine read_rest_of_file(z,dz,rhoc,nxc,nyc,nzc,FileName)
  
      integer :: coordsys
      integer :: nxc,nyc,nzc
      integer :: nvar,nscalars
      integer :: selfgrav_boolean, particles_boolean
      real(c_double)       :: gamma1c,csc
      real(c_double)       :: tc,dtc
      real(c_double), dimension(nxc) :: xc
      real(c_double), dimension(nyc) :: yc
      real(c_double), dimension(nzc) :: zc
      real(c_double), dimension(nzc,nyc,nxc) :: ruxc,ruyc,ruzc,engc
      character(len=90)   :: FileName
      integer :: i
  
      real, intent(out) :: dz
      real, dimension(mz), intent(out) :: z
      real, dimension(nz,ny,nx), intent(out) :: rhoc
  
  !order: 
  !coordsys,nx,ny,nz,nvar,nscalars,selfgrav_boolean, particles_boolean,gamma1,cs,t,dt,x,y,z,rho,rux,ruy,ruz,eng
  !
  
      open(99, file = trim(FileName), form = 'unformatted', ACCESS='stream')
      read(99) coordsys
      read(99) nxc,nyc,nzc,nvar,nscalars
      read(99) selfgrav_boolean, particles_boolean
      read(99) gamma1c,csc
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
    endsubroutine read_rest_of_file
!************************************************************************************
endmodule ReadC
