program flux_feautrier

  use Auxiliary
  use Common
  use Columns
  use Disk
  use GasState
  use ContinuousOpacity
  use ReadAthena
  use FileIO
  
  implicit none

  type (column_case) :: c
  
  real, dimension(mz,nyloc,nxloc,nw) :: U
  real, dimension(nz,nyloc,nxloc,nw) :: absorp_coeff
  real, dimension(nz,nyloc,nxloc) :: rho3d,temp3d
  !real, dimension(nz,nw) :: V,Ip,Im
!
  real, dimension(mz) :: z
  real, dimension(nz) :: aa,bb,cc,dd
  real, dimension(nz) :: opacity, albedo
  real, dimension(nz) :: source_function

  real, dimension(nw) :: waves_angstrom
!
  real :: dz,z0,z1
  real :: start, finish
  real :: w0=3000,w1=5000
  real :: sigma_grey
!
  integer :: iw,ix,iy,iprocx,iprocy
!
  logical :: lgrey=.false.,lread_athena=.true.
  character(len=90) :: snapshot
  character(len=90) :: inputfile='./input.in'
!
  !integer, parameter :: ncolumns=4

  !character (len=90), parameter, dimension(ncolumns) :: column_names = &
  !(/'rho','T','rho1','T1'/)
!
  namelist /input/ z0,z1,w0,w1,sigma_grey,lgrey,lread_athena
!
! Start the time counter
!
  call cpu_time(start)
!
  !c%T = 0.
  !c%rho = 0.
  !c%rho1 = 0.
  !c%T1 = 0.
!
!  Read the input namelist with the user-defined parameters.
!
  call getarg(1,snapshot)
  if (snapshot=='') snapshot="0000"      
!
  open(20,file=trim(inputfile))
  read(20,nml=input)
  close(20)
!
! Do a sanity check for grey RT
!
  if (lgrey.and.(nw /= 1)) then
    print*,"For Grey RT use only one wavelength. Switch nw=1 in resolution.in"
    stop
  endif
!
! Calculate the grid variables
!
  call read_temperature_input(inputfile)
  call read_density_input(inputfile)
  call read_gas_state_input(inputfile)
!
  call calc_wavelength(w1,w0,waves_angstrom)
!
  call pre_calc_opacity_quantities(waves_angstrom)
! 
  if (lread_athena) call read_athena_input(inputfile)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through processors  
  do iprocx=0,nprocx-1; do iprocy=0,nprocy-1   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  lroot = (iprocx==0) .and. (iprocy==0)
  if (lread_athena) then
    call read_from_athena(z,dz,rho3d,temp3d,iprocx,iprocy,snapshot)
  else
    call calc_grid(z1,z0,z,dz)
    call calc_density(c%rho,z)
    call calc_temperature(c%T,z)
  endif
!
!***********************************************************************
!xy-dependent starts here
!***********************************************************************
  xloop: do ix=1,nxloc
    yloop: do iy=1,nyloc
      if (lread_athena) then
        c%rho =  rho3d(1:nz,iy,ix)
        c%T   = temp3d(1:nz,iy,ix)
      endif
      c%rho1=1./c%rho
      c%T1=1./c%T
      call calc_hydrogen_ion_frac(c)
      call solve_gas_state(c)!rho,rho1,NHII_NHINHII,number_density,inv_number_density,nHI,nHII,ne,ionization_factor)
      c%electron_pressure = get_electron_pressure(c) !ne,T)
      c%e_scatter         = get_electron_thomson_scattering(c)!inv_number_density,nHII)
      c%theta             = 5040.* c%T1
      c%theta1            = 1./c%theta
      c%lgtheta           = log10(c%theta)
      c%lgtheta2          = c%lgtheta**2
      c%hm_bf_factor      = get_hydrogen_ion_bound_free(c) !electron_pressure,theta)
!
! Loop over wavelengths
!
      wavelength: do iw=1,nw
        lfirst=lroot.and.(ix==1).and.(iy==1).and.(iw==1) 
!    
        if (lgrey) then
          call grey_parameters(c,sigma_grey,source_function,opacity,albedo) !rho,T,sigma_grey,source_function,opacity,albedo)
        else
          source_function = get_source_function(c%T1,iw)
          c%stim_factor = get_hydrogen_stimulated_emission(iw,c%theta)
          call calc_opacity_and_albedo(c,iw,opacity,albedo)
        endif
!
! Populate coefficient arrays
!
        call fill_center_coeffs(aa,bb,cc,dd,opacity,albedo,source_function,dz)
        call fill_boundary_coeffs(aa,bb,cc,dd,opacity,albedo,source_function,dz)    
!
! Solve the system of equations
!
        if (lfirst) print*,'sum(a), sum(b), sum(c), sum(d)=',sum(aa), sum(bb), sum(cc), sum(dd)
        call tridag(aa, bb, cc, dd, U(n1:n2,iy,ix,iw))
        call update_ghosts(U(:,iy,ix,iw))
        if (lfirst) print*, 'min/max(U)', minval(U(n1:n2,iy,ix,iw)), maxval(U(n1:n2,iy,ix,iw))
!
        absorp_coeff(:,iy,ix,iw)=opacity
!
      enddo wavelength
    enddo yloop
  enddo xloop
!
! Calculate post-processing on one column for diagnostic: flux and intensities.
!
  !call calc_auxiliaries(U(:,ny,nx,:),absorp_coeff(:,ny,nx,:),dz,V,Ip,Im)
!
! Write output
!
  if (lroot) then
    call output_grid(z)
    !Output for diagnostic purposes 
    call output_ascii(U(:,nyloc,nxloc,:),absorp_coeff(:,nyloc,nxloc,:))
  endif
  call output_binary(U,absorp_coeff,iprocx,iprocy,snapshot)
  enddo;enddo
! 
! Finish the time counter and print the wall time
!
  call cpu_time(finish)
  print*,"Wall time = ",finish-start," seconds."
  print*,"Execution time =",(finish-start)/(nz*ny*nx*nw)*1e6,&
       " micro-seconds per frequency point per mesh point." 
!
endprogram flux_feautrier
!***********************************************************************
