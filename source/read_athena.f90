module ReadAthena

  use Common
  use Disk
  use fileIO
  use, intrinsic :: iso_c_binding, only: c_float, c_double, c_int  
  implicit none

  private

  public :: read_from_athena

  character(len=90)   :: RunName
  character(len=90)   :: datadir="./output_athena/distributed"
  character(len=90)   :: snapshot="0000"

  real :: Mbh_SolarMasses,r0ref_rg
  real :: aspect_ratio,mean_molecular_weight,rho0

  integer :: nprocz=2,nprocy=1,nprocx=2
  integer :: iprocx=0,iprocy=0
  
  namelist /athena_input/ RunName,Mbh_SolarMasses,r0ref_rg,&
       aspect_ratio,mean_molecular_weight,rho0,datadir,snapshot,&
       nprocz,nprocy,nprocx,iprocx,iprocy

  
contains
!************************************************************************************
    subroutine read_from_athena(z,dz,rho,temp)
!
      integer(c_int) :: coordsys
      integer(c_int) :: nxloc_,nyloc_,nzloc
      integer(c_int) :: nvar,nscalars
      integer(c_int) :: selfgrav_boolean, particles_boolean
      real(c_double) :: gamma1,cs
      real(c_double) :: t,dt
!      
      real(c_double), allocatable :: xloc(:),yloc(:),zloc(:)
      real(c_double), allocatable :: tmp_loc(:,:,:)
!
      real, dimension(nz,nyloc,nxloc), intent(out) :: rho,temp
!
      real, dimension(mz), intent(out) :: z
      real, dimension(nz) :: zn      
      real, intent(out) :: dz
      integer :: iproc,iprocz
      character(len=90)   :: head,tail,filename,sproc,base
!
      real :: z0
      integer :: iz0
      integer :: iz1

  !order: 
  !coordsys,nx,ny,nz,nvar,nscalars,selfgrav_boolean, particles_boolean,gamma1,cs,t,dt,x,y,z,rho,rux,ruy,ruz,eng
  !
      open(40,file='./input.in')
      read(40,nml=athena_input)
      close(40)

      base=trim(datadir)//"/id"
      tail="."//trim(snapshot)//".bin"

      do iprocz=0,nprocz-1

        iproc = modulo(iprocz,nprocz) * nprocx*nprocy + modulo(iprocy,nprocy) * nprocx + modulo(iprocx,nprocx)

        sproc=itoa(iproc)
        head=trim(base)//trim(sproc)//"/"//trim(RunName)
        if (iproc==0) then 
          filename= trim(head)//trim(tail)
        else
          filename= trim(head)//"-id"//trim(sproc)//trim(tail)
        endif
        if (file_exists(trim(filename)) .eqv. .false.) then
          print*,"iproc=",iproc
          print*,"File ",trim(filename)," does not exist"
          stop
        endif

        print*,"Reading ",trim(FileName)
        open(99, file = trim(FileName), form = 'unformatted', ACCESS='stream')

        read(99) coordsys
        read(99) nxloc_,nyloc_,nzloc,nvar,nscalars
        read(99) selfgrav_boolean, particles_boolean
        read(99) gamma1,cs
        read(99) t,dt

        allocate(xloc(nxloc_)); read(99) xloc
        allocate(yloc(nyloc_)); read(99) yloc
        allocate(zloc(nzloc)); read(99) zloc
!
        if (iprocz==0) then
            z0=zloc(1)
            dz=zloc(2)-zloc(1)
         endif

        iz0 = nint((zloc(1)-z0)/dz) + 1
        iz1 = iz0 + nzloc - 1
!
        zn(iz0:iz1) = zloc
!
! Sanity check
!
        if (nxloc_ /= nxloc) then
          print*,"Resolution in x from Athena = ",nxloc_
          print*,"RT x-dimensionality from resolution.in = ",nxloc
          print*,"The Athena x dimensionality does not match the chosen for the RT post-processing. Fix."
          stop
        endif
        if (nyloc_ /= nyloc) then
          print*,"Resolution in y from Athena = ",nyloc_
          print*,"RT y-dimensionality from resolution.in = ",nyloc
          print*,"The Athena y dimensionality does not match the chosen for the RT post-processing. Fix."
          stop
        endif
        if (iprocz==nprocz-1) then
          if (iz1 /= nz) then
            print*,"Resolution in z from Athena = ",iz1
            print*,"RT z-dimensionality from resolution.in = ",nz
            print*,"The Athena z dimensionality does not match the chosen for the RT post-processing. Fix."
            stop
          endif
        endif
!
        allocate(tmp_loc(nzloc,nyloc_,nxloc_))
!
        !rho1=1./rho(:,iy,ix)
        !ekin = 0.5*rho1*ru2(:,iy,ix)
        !eint = eng(:,iy,ix) - ekin
        !cs2 = gamma*gamma1 * eint * rho1
        !temp(:,iy,ix) = mean_molecular_weight * amu * cs2 * gamma_inv * k1_cgs

        read(99) tmp_loc
        rho(iz0:iz1,:,:) = tmp_loc
!
! Build kinetic energy; ekin = (rux**2+ruy**2+ruz**2)/(2*rho)
!        
        read(99) tmp_loc !rux
        temp(iz0:iz1,:,:) = .5*tmp_loc**2/rho(iz0:iz1,:,:)
        read(99) tmp_loc !ruy
        temp(iz0:iz1,:,:) = temp(iz0:iz1,:,:) + .5*tmp_loc**2/rho(iz0:iz1,:,:)
        read(99) tmp_loc !ruz
        temp(iz0:iz1,:,:) = temp(iz0:iz1,:,:) + .5*tmp_loc**2/rho(iz0:iz1,:,:)
!        
! Define internal energy; eint = eng - ekin.
! Then the sound speed; cs2 = gamma*gamma1 * eint * rho1
! The gammas and the constants will be added in the next subroutine 
!        
        read(99) tmp_loc !eng
        temp(iz0:iz1,:,:) = (tmp_loc-temp(iz0:iz1,:,:))/rho(iz0:iz1,:,:)
        
        deallocate(tmp_loc,xloc,yloc,zloc)
        
        close(99)        
      enddo
!  
      print*,"Done reading all files"
!
      call postprocess_athena_values(rho,temp,zn,z,dz,gamma1)
!
    endsubroutine read_from_athena
!************************************************************************************
    subroutine postprocess_athena_values(rho,temp,zn,z,dz,gamma1)
!
      real, dimension(nz,nyloc,nxloc), intent(inout) :: rho, temp
      real, dimension(mz), intent(out) :: z
      real, dimension(nz) :: zn
      real, intent(in) :: gamma1
      real, intent(out) :: dz
!
      real :: rg,Mbh,rr,g0,Omega,H
      real :: unit_time,unit_length,unit_velocity
      real :: unit_density,unit_temperature
      integer :: i,ix,iy
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
      unit_density  = rho0
      unit_time     = 1./Omega
      unit_length   = H
      unit_velocity = unit_length/unit_time
      unit_temperature = gamma1 * mean_molecular_weight * amu *  k1_cgs  * unit_velocity**2      
!
      rho = rho*unit_density      
      temp = temp*unit_temperature
!
      zn=zn*unit_length
      dz=(zn(nz)-zn(1))/(nz-1)
      z(n1:n2)=zn
! Fill in ghost zones
      do i=1,ng
        z(n1-i) = z(n1) - i*dz
        z(n2+i) = z(n2) + i*dz
      enddo
!
      do ix=1,nxloc; do iy=1,nyloc
        call calc_temperature(temp(:,iy,ix),z,lfrom_read_athena=.true.)
      enddo; enddo
!
      print*,'unit_time        = ',unit_time,' s'
      print*,'unit_length      = ',unit_length,' cm'
      print*,'unit_velocity    = ',unit_velocity,' cm/s'
      print*,'unit_density     = ',unit_density,' g/cm3'
      print*,'unit_temperature = ',unit_temperature,' K'
!
      print*,'minval(rho),maxval(rho),mean(rho)',&
           minval(rho),maxval(rho),sum(rho)/(nxloc*nyloc*nz)
      print*,'minval(T),maxval(T),mean(T)',&
           minval(temp),maxval(temp),sum(temp)/(nxloc*nyloc*nz)
!
    endsubroutine postprocess_athena_values
!************************************************************************************
    character (len=21) function itoa(n)
! plucked from pencil
      integer, intent(in) :: n
!
      write (itoa, '(I21)') n   ! (64 bit integer plus sign)
      itoa = adjustl(itoa)
!
    endfunction itoa
!************************************************************************************
  endmodule ReadAthena
