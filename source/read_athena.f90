module ReadAthena

  use Common
  use Disk
  use FileIO
  use, intrinsic :: iso_c_binding, only: c_float, c_double, c_int  
  implicit none

  private

  public :: read_from_athena,read_athena_input

  character(len=90)   :: RunName
  character(len=90)   :: datadir="./input_athena/distributed"

  real :: Mbh_SolarMasses,r0ref_rg
  real :: aspect_ratio,mean_molecular_weight,rho0,rho_floor

  real :: rg,Mbh,rr,g0,Omega,H
  real :: unit_time,unit_length,unit_velocity
  real :: unit_density,unit_temperature

    namelist /athena_input/ RunName,Mbh_SolarMasses,r0ref_rg,&
      aspect_ratio,mean_molecular_weight,rho0,datadir,rho_floor
  
contains
!************************************************************************************
  subroutine read_athena_input(inputfile)
!
    character(len=90)   :: inputfile
!    
    open(40,file=trim(inputfile))
    read(40,nml=athena_input)
    close(40)
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
    unit_temperature = mean_molecular_weight * amu *  k1_cgs  * unit_velocity**2
!
    print*,'unit_time        = ',unit_time,' s'
    print*,'unit_length      = ',unit_length,' cm'
    print*,'unit_velocity    = ',unit_velocity,' cm/s'
    print*,'unit_density     = ',unit_density,' g/cm3'
    print*,'unit_temperature = ',unit_temperature,' K'
!
  endsubroutine read_athena_input
!************************************************************************************
  subroutine read_from_athena(z,dz,rho,temp,iprocx,iprocy,snapshot)
!
      integer(c_int) :: nxloc_,nyloc_,nzloc
      real(c_double) :: gamma1

      integer(c_int) :: dummy1_int ! coordsys
      integer(c_int) :: dummy2_int, dummy3_int !nvar,nscalars
      integer(c_int) :: dummy4_int, dummy5_int !selfgrav_boolean, particles_boolean
      real(c_double) :: dummy1_real,dummy2_real,dummy3_real !cs,t,dt
!      
      real(c_double), dimension(nxloc) :: xloc
      real(c_double), dimension(nyloc) :: yloc
      real(c_double), allocatable :: zloc(:)
      !real(c_double), dimension(nzloc,nyloc,nxloc) :: tmp_loc

!      real(c_double), allocatable :: xloc(:),yloc(:),zloc(:)
      real(c_double), allocatable :: tmp_loc(:,:,:),tmp_loc_transpose(:,:,:)
!
      real, dimension(nz,nyloc,nxloc), intent(out) :: rho,temp
!
      real, dimension(mz), intent(out) :: z
      real, dimension(nz) :: zn      
      real, intent(out) :: dz
      integer :: iproc,iprocz
      integer, intent(in) :: iprocx,iprocy
      character(len=90)   :: head,tail,filename,sproc,base,snapshot
!
      real :: z0
      integer :: iz0
      integer :: iz1

  !order: 
  !coordsys,nx,ny,nz,nvar,nscalars,selfgrav_boolean, particles_boolean,gamma1,cs,t,dt,x,y,z,rho,rux,ruy,ruz,eng
  !
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

        open(99, file = trim(filename), form = 'unformatted', ACCESS='stream')

        read(99) dummy1_int !coordsys
!
! Sanity check
!
        read(99) nxloc_
        if (nxloc_ /= nxloc) then
          print*,"Resolution in x from Athena = ",nxloc_
          print*,"RT x-dimensionality from resolution.in = ",nxloc
          print*,"The Athena x dimensionality does not match the chosen for the RT post-processing. Fix."
          stop
        endif
        read(99) nyloc_
        if (nyloc_ /= nyloc) then
          print*,"Resolution in y from Athena = ",nyloc_
          print*,"RT y-dimensionality from resolution.in = ",nyloc
          print*,"The Athena y dimensionality does not match the chosen for the RT post-processing. Fix."
          stop
        endif
        read(99) nzloc
!       
        read(99) dummy2_int,dummy3_int,dummy4_int,dummy5_int !nvar,nscalars,selfgrav_boolean, particles_boolean
        read(99) gamma1
        read(99) dummy1_real,dummy2_real,dummy3_real !cs,t,dt
!        
        read(99) xloc
        read(99) yloc
        allocate(zloc(nzloc))
        read(99) zloc
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
        allocate(tmp_loc(nxloc,nyloc,nzloc))
        allocate(tmp_loc_transpose(nzloc,nyloc,nxloc))
!
        !rho1=1./rho(:,iy,ix)
        !ekin = 0.5*rho1*ru2(:,iy,ix)
        !eint = eng(:,iy,ix) - ekin
        !cs2 = gamma*gamma1 * eint * rho1
        !temp(:,iy,ix) = mean_molecular_weight * amu * cs2 * gamma_inv * k1_cgs

        read(99) tmp_loc; call transpose_loc(tmp_loc,tmp_loc_transpose)
        rho(iz0:iz1,:,:) = tmp_loc_transpose
!
! Build kinetic energy; ekin = (rux**2+ruy**2+ruz**2)/(2*rho)
!        
        read(99) tmp_loc; call transpose_loc(tmp_loc,tmp_loc_transpose) !rux
        temp(iz0:iz1,:,:) = .5*tmp_loc_transpose**2/rho(iz0:iz1,:,:)
        read(99) tmp_loc; call transpose_loc(tmp_loc,tmp_loc_transpose) !ruy
        temp(iz0:iz1,:,:) = temp(iz0:iz1,:,:) + .5*tmp_loc_transpose**2/rho(iz0:iz1,:,:)
        read(99) tmp_loc; call transpose_loc(tmp_loc,tmp_loc_transpose) !ruz
        temp(iz0:iz1,:,:) = temp(iz0:iz1,:,:) + .5*tmp_loc_transpose**2/rho(iz0:iz1,:,:)
!        
! Define internal energy; eint = eng - ekin.
! Then the sound speed; cs2 = gamma*gamma1 * eint * rho1
! The gammas and the constants will be added in the next subroutine 
!        
        read(99) tmp_loc; call transpose_loc(tmp_loc,tmp_loc_transpose) !eng
        temp(iz0:iz1,:,:) = (tmp_loc_transpose-temp(iz0:iz1,:,:))/rho(iz0:iz1,:,:)
        
        deallocate(tmp_loc,zloc,tmp_loc_transpose)
        
        close(99)        
      enddo
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
      integer :: i,ix,iy
!
! Apply density floor after converting to physical units
!
      rho = rho*unit_density
      if (ISNAN(rho_floor)) then
        print*,""
        print*,"Must set `rho_floor` (cgs) value in namelist: `input.in/athena_input`."
        print*,""
        stop ! quit program
      else
        ! Apply density floor
        where (rho < rho_floor) rho = rho_floor
      endif
!
! convert temperature to physical units
!
      temp = temp*gamma1*unit_temperature
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
      if (lroot) then
        print*,'minval(rho),maxval(rho),mean(rho)',&
             minval(rho),maxval(rho),sum(rho)/(nxloc*nyloc*nz)
        print*,'minval(T),maxval(T),mean(T)',&
             minval(temp),maxval(temp),sum(temp)/(nxloc*nyloc*nz)
      endif
!
    endsubroutine postprocess_athena_values
!************************************************************************************
    subroutine transpose_loc(a,b)
!
      real, dimension(:,:,:),intent(in) :: a
      real, dimension(:,:,:),intent(out) :: b
      real, dimension(size(a,1),size(a,2)) :: buffer
      integer :: i
!
      do i=1,size(a,3)
         buffer   = a(:,:,i)
         b(i,:,:) = transpose(buffer)
      enddo
!
    endsubroutine transpose_loc
!************************************************************************************
  endmodule ReadAthena
