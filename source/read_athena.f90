module ReadAthena

  use Common
  use Disk
  use fileIO
  use, intrinsic :: iso_c_binding, only: c_float, c_double, c_int  
  implicit none

  private

  public :: read_from_athena

  character(len=90)   :: RunName
  character(len=90)   :: FileName
  
  real :: Mbh_SolarMasses,r0ref_rg
  real :: aspect_ratio,mean_molecular_weight,rho0

  integer :: nproc=1
  
  namelist /athena_input/ RunName,Mbh_SolarMasses,r0ref_rg,&
       aspect_ratio,mean_molecular_weight,rho0,nproc

  
contains
!************************************************************************************
    subroutine read_from_athena(z,dz,rho,temp)
  
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
      !!!!!
      
      real, intent(out) :: dz
      real, dimension(mz), intent(out) :: z
      real, dimension(nz,ny,nx), intent(out) :: rho,temp
      real, dimension(nz,ny,nx) :: rux,ruy,ruz,eng
      real, dimension(nz) :: zn      
      integer :: iproc
      character(len=90)   :: Head,Tail

  !order: 
  !coordsys,nx,ny,nz,nvar,nscalars,selfgrav_boolean, particles_boolean,gamma1,cs,t,dt,x,y,z,rho,rux,ruy,ruz,eng
  !
      open(40,file='./input.in')
      read(40,nml=athena_input)
      close(40)

      Head="./output_athena/distributed/id"
      Tail=".0000.bin"
      
      do iproc=0,nproc-1
         
        if (iproc==0) then 
          FileName= trim(Head)//trim(itoa(iproc))//"/"//trim(RunName)//trim(Tail)
        else
           FileName= trim(Head)//trim(itoa(iproc))//"/"//trim(RunName)//"-id"//trim(itoa(iproc))//trim(Tail)
        endif
        if (file_exists(trim(Filename)) .eqv. .false.) then
           print*,"File ",trim(Filename)," does not exist"
           stop
        endif

        print*,"iproc=",iproc
        print*,trim(FileName)
        print*,''        
        open(99, file = trim(FileName), form = 'unformatted', ACCESS='stream')

        
        read(99) coordsys
        print*,'coordsys',coordsys
        print*,''
        
        read(99) nxloc,nyloc,nzloc,nvar,nscalars
        print*,'nxloc,nyloc,nzloc,nvar,nscalars',nxloc,nyloc,nzloc,nvar,nscalars
        print*,''
        
        read(99) selfgrav_boolean, particles_boolean
        print*,'selfgrav_boolean, particles_boolean',selfgrav_boolean, particles_boolean
        print*,''
        
        read(99) gamma1,cs
        print*,'gamma1,cs',gamma1,cs
        print*,''
        
        read(99) t,dt
        print*,'t, dt',t, dt        
        print*,''

        allocate(xloc(nxloc))
        read(99) xloc
        print*,'xloc=',xloc
        print*,''
        deallocate(xloc)

        allocate(yloc(nyloc))
        read(99) yloc
        print*,'yloc=',yloc
        print*,''
        deallocate(yloc)

        allocate(zloc(nzloc))        
        read(99) zloc
        print*,'zloc=',zloc
        print*,''
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
        print*,'rho=',rho_loc(:,1,1)
        print*,''
        if (nproc==1) rho=rho_loc
        deallocate(rho_loc)

        allocate(rux_loc(nzloc,nyloc,nxloc))
        read(99) rux_loc
        print*,'rux_loc=',rux_loc(:,1,1)
        print*,''
        if (nproc==1) rux=rux_loc
        deallocate(rux_loc)

        allocate(ruy_loc(nzloc,nyloc,nxloc))
        read(99) ruy_loc
        print*,'ruy_loc=',ruy_loc(:,1,1)
        print*,''
        if (nproc==1) ruy=ruy_loc
        deallocate(ruy_loc)
        
        allocate(ruz_loc(nzloc,nyloc,nxloc))
        read(99) ruz_loc
        print*,'ruz_loc=',ruz_loc(:,1,1)
        print*,''
        if (nproc==1) ruz=ruz_loc
        deallocate(ruz_loc)
        
        allocate(eng_loc(nzloc,nyloc,nxloc))
        read(99) eng_loc
        print*,'eng_loc=',eng_loc(:,1,1)        
        print*,''
        if (nproc==1) eng=eng_loc
        deallocate(eng_loc)
        
        print*,'!!!!!!!!!!!!!!!!'

        print*,''
        
        close(99)        
        !print*,'closed'
      enddo
!  
      print*,'done reading'
      !stop
!
      call postprocess_athena_values(zn,gamma1,rho,rux,ruy,ruz,eng,z,temp,dz)
!
    endsubroutine read_from_athena
!************************************************************************************
    subroutine postprocess_athena_values(zn,gamma1,rho,rux,ruy,ruz,eng,z,temp,dz)

      integer :: i,ix,iy
      integer :: nxc,nyc,nzc
      real :: gamma,gamma_inv,rg,Mbh,rr,g0,Omega,H
      real :: unit_time,unit_length,unit_velocity
      real :: unit_density,unit_mass,unit_energy,unit_edens
      real, intent(in) :: gamma1
      real, intent(out) :: dz
      real, dimension(nz) :: zn
      real, dimension(nz) :: ekin,eint,cs2,rho1
      real, dimension(mz), intent(out) :: z
      real, dimension(nz,ny,nx) :: rux,ruy,ruz,eng
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
!
! Fill in ghost zones
!
      do i=1,ng
         z(n1-i) = z(n1) - i*dz
         z(n2+i) = z(n2) + i*dz
      enddo
!
      rho = rho*unit_density
!
      rux = rux*unit_density*unit_velocity
      ruy = ruy*unit_density*unit_velocity
      ruz = ruz*unit_density*unit_velocity
!
      gamma = 1+gamma1
      gamma_inv = 1./gamma
!
      eng = eng*unit_energy/unit_length**3
      do ix=1,nx
        do iy=1,ny
          rho1=1./rho(:,iy,ix)
          ekin = 0.5*rho1*(rux(1:nz,iy,ix)**2 + &
                           ruy(1:nz,iy,ix)**2 + &
                           ruz(1:nz,iy,ix)**2)
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
