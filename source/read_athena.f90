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
  
  real :: Mbh_SolarMasses,r0ref_rg
  real :: aspect_ratio,mean_molecular_weight,rho0

  integer :: nproc=1
  
  namelist /athena_input/ RunName,Mbh_SolarMasses,r0ref_rg,&
       aspect_ratio,mean_molecular_weight,rho0,nproc,datadir

  
contains
!************************************************************************************
    subroutine read_from_athena(z,dz,rho,temp)
!
      integer(c_int) :: coordsys
      integer(c_int) :: nxloc,nyloc,nzloc
      integer(c_int) :: nvar,nscalars
      integer(c_int) :: selfgrav_boolean, particles_boolean
      real(c_double)       :: gamma1,cs
      real(c_double)       :: t,dt
      real(c_double), allocatable :: xloc(:)
      real(c_double), allocatable :: yloc(:)
      real(c_double), allocatable :: zloc(:)
      real(c_double), allocatable :: rho_loc(:,:,:),rux_loc(:,:,:)
      real(c_double), allocatable :: ruy_loc(:,:,:),ruz_loc(:,:,:),eng_loc(:,:,:)
!
      real, intent(out) :: dz
      real, dimension(mz), intent(out) :: z
      real, dimension(nz,ny,nx), intent(out) :: rho,temp
      real, dimension(nz,ny,nx) :: ru2,eng
      real, dimension(nz) :: zn      
      integer :: iproc
      character(len=90)   :: head,tail,filename,sproc,base

  !order: 
  !coordsys,nx,ny,nz,nvar,nscalars,selfgrav_boolean, particles_boolean,gamma1,cs,t,dt,x,y,z,rho,rux,ruy,ruz,eng
  !
      open(40,file='./input.in')
      read(40,nml=athena_input)
      close(40)

      base=trim(datadir)//"/id"
      tail=".0000.bin"
      
      do iproc=0,nproc-1
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
        read(99) nxloc,nyloc,nzloc,nvar,nscalars
        read(99) selfgrav_boolean, particles_boolean
        read(99) gamma1,cs
        read(99) t,dt
        
        allocate(xloc(nxloc)); read(99) xloc; deallocate(xloc)
        allocate(yloc(nyloc)); read(99) yloc; deallocate(yloc)

        allocate(zloc(nzloc))        
        read(99) zloc
        if (nproc==1) zn=zloc
        deallocate(zloc)
        
        !print*,"Athena input x dimensionality=",nxc
        !if (nxc /= nx) then
        !   print*,"RT x-dimensionality from resolution.in =",nx
        !   print*,"The Athena x dimensionality does not match the chosen for the RT post-processing. Fix."
        !   print*,''           
        !   !stop
        !endif
!
        !print*,"Athena input y dimensionality=",nyc
        !if (nyc /= ny) then
        !   print*,"RT y-dimensionality from resolution.in =",ny
        !   print*,"The Athena y dimensionality does not match the chosen for the RT post-processing. Fix."
        !   print*,''           
        !   !stop
        !endif
!
        !print*,"Athena input z dimensionality=",nzc
        !if (nzc /= nz) then
        !   print*,"RT z-dimensionality from resolution.in =",nz
        !   print*,"The Athena z dimensionality does not match the chosen for the RT post-processing. Fix."
        !   print*,''           
        !   !stop
        !endif

        allocate(rho_loc(nzloc,nyloc,nxloc))
        read(99) rho_loc
        if (nproc==1) rho=rho_loc
        deallocate(rho_loc)

        allocate(rux_loc(nzloc,nyloc,nxloc))
        read(99) rux_loc
        if (nproc==1) ru2=rux_loc**2
        deallocate(rux_loc)

        allocate(ruy_loc(nzloc,nyloc,nxloc))
        read(99) ruy_loc
        if (nproc==1) ru2=ru2+ruy_loc**2
        deallocate(ruy_loc)
        
        allocate(ruz_loc(nzloc,nyloc,nxloc))
        read(99) ruz_loc
        if (nproc==1) ru2=ru2+ruz_loc**2
        deallocate(ruz_loc)
        
        allocate(eng_loc(nzloc,nyloc,nxloc))
        read(99) eng_loc
        if (nproc==1) eng=eng_loc
        deallocate(eng_loc)
        
        close(99)        
      enddo
!  
      print*,"Done reading all files"
!
      call postprocess_athena_values(zn,gamma1,rho,ru2,eng,z,temp,dz)
!
    endsubroutine read_from_athena
!************************************************************************************
    subroutine postprocess_athena_values(zn,gamma1,rho,ru2,eng,z,temp,dz)

      integer :: i,ix,iy
      real :: gamma,gamma_inv,rg,Mbh,rr,g0,Omega,H
      real :: unit_time,unit_length,unit_velocity
      real :: unit_density,unit_mass,unit_energy,unit_edens
      real, intent(in) :: gamma1
      real, intent(out) :: dz
      real, dimension(nz) :: zn
      real, dimension(nz) :: ekin,eint,cs2,rho1
      real, dimension(mz), intent(out) :: z
      real, dimension(nz,ny,nx) :: ru2,eng
      real, dimension(nz,ny,nx), intent(inout) :: rho
      real, dimension(nz,ny,nx), intent(out) :: temp
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
      zn=zn*unit_length
!      
      dz=(zn(nz)-zn(1))/(nz-1)
      z(n1:n2)=zn
! Fill in ghost zones
      do i=1,ng
         z(n1-i) = z(n1) - i*dz
         z(n2+i) = z(n2) + i*dz
      enddo
!
      rho = rho*unit_density
      ru2 = ru2*(unit_density*unit_velocity)**2
      eng = eng*unit_energy/unit_length**3      
!
      gamma = 1+gamma1
      gamma_inv = 1./gamma
!
      do ix=1,nx
        do iy=1,ny
          rho1=1./rho(:,iy,ix)
          ekin = 0.5*rho1*ru2(:,iy,ix)
          eint = eng(:,iy,ix) - ekin
          cs2 = gamma*gamma1 * eint * rho1
          temp(:,iy,ix) = mean_molecular_weight * amu * cs2 * gamma_inv * k1_cgs
!
          call calc_temperature(temp(:,iy,ix),z,&
               lfrom_read_athena=.true.)
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
           minval(rho),maxval(rho),sum(rho)/(nx*ny*nz)
      print*,'minval(T),maxval(T),mean(T)',&
           minval(temp),maxval(temp),sum(temp)/(nx*ny*nz)
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
