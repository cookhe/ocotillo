program flux_feautrier

  use Auxiliary
  use Common
  use Disk
  use GasState
  use ContinuousOpacity
  use ReadAthena
  use FileIO
  
  implicit none

  real, dimension(mz,nyloc,nxloc,nw) :: U
  real, dimension(nz,nyloc,nxloc,nw) :: absorp_coeff
  real, dimension(nz,nyloc,nxloc) :: rho3d,temp3d
  !real, dimension(nz,nw) :: V,Ip,Im
!
  real, dimension(mz) :: z
  real, dimension(nz) :: aa,bb,cc,dd
  real, dimension(nz) :: rho,rho1,T,T1
  real, dimension(nz) :: opacity, albedo
  real, dimension(nz) :: NHII_NHINHII,number_density,inv_number_density
  real, dimension(nz) :: nHII,nHI,ne,ionization_factor
  real, dimension(nz) :: e_scatter,theta,theta1,lgtheta,lgtheta2
  real, dimension(nz) :: electron_pressure,hm_bf_factor
  real, dimension(nz) :: stim_factor,source_function

  real, dimension(nw) :: waves_cm,waves1_cm,nu_Hz
  real, dimension(nw) :: waves_angstrom,waves1_angstrom
!
  real :: wave_cm,wave_angstrom,wave1_cm,wave1_angstrom
  real :: dz,z0,z1
  real :: start, finish
  real :: w0=3000,w1=5000
  real :: sigma_grey
!
  integer :: iw,ix,iy,iprocx,iprocy
  integer :: log_overflow_limit
!
  logical :: lgrey=.false.,lread_athena=.true.
  logical :: lroot
  character(len=90) :: snapshot
  character(len=90) :: inputfile='./input.in'
!
  namelist /input/ z0,z1,w0,w1,sigma_grey,lgrey,lread_athena
!
! Start the time counter
!
  call cpu_time(start)
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
  log_overflow_limit = int(floor(log10(float_info_max)))
!
! Calculate the grid variables
!
  call read_temperature_input(inputfile)
  call read_density_input(inputfile)
  call read_gas_state_input(inputfile)
!
  call calc_wavelength(w1,w0,waves_angstrom,&
       waves1_angstrom,waves_cm,waves1_cm,nu_Hz)
!
  call pre_calc_opacity_quantities()
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
     call calc_density(rho,z)
     call calc_temperature(T,z)
  endif
!
!***********************************************************************
!xy-dependent starts here
!***********************************************************************
  xloop: do ix=1,nxloc
    yloop: do iy=1,nyloc
      if (lread_athena) then
        rho = rho3d(1:nz,iy,ix)
        T   = temp3d(1:nz,iy,ix)
      endif
      rho1=1./rho
      T1=1./T
      call calc_hydrogen_ion_frac(rho1,T,T1,NHII_NHINHII)
      call solve_gas_state(rho,rho1,NHII_NHINHII,number_density,inv_number_density,nHI,nHII,ne,ionization_factor)
      electron_pressure = get_electron_pressure(ne,T)
      e_scatter         = get_electron_thomson_scattering(inv_number_density,nHII)
      theta             = 5040.* T1
      theta1            = 1./theta
      lgtheta           = log10(theta)
      lgtheta2          = lgtheta**2
      hm_bf_factor      = get_hydrogen_ion_bound_free(electron_pressure,theta)
!
! Loop over wavelengths
!
      wavelength: do iw=1,nw
         lfirst=(ix==1).and.(iy==1).and.(iw==1) 
         wave_cm = waves_cm(iw)
         wave1_cm = waves1_cm(iw)
         if (lfirst) print*, 'waves_cm min/max',minval(waves_cm),maxval(waves_cm)
         wave_angstrom = waves_angstrom(iw)
         wave1_angstrom = waves1_angstrom(iw)
!    
         if (lgrey) then
            call grey_parameters(rho,T,sigma_grey,source_function,opacity,albedo)
         else
            source_function = get_source_function(wave1_cm,T1,log_overflow_limit)
            stim_factor = get_hydrogen_stimulated_emission(wave1_angstrom,theta)
            call calc_opacity_and_albedo(e_scatter,rho,rho1,ne,NHII_NHINHII,nHI,nHII,&
                 T,T1,theta,theta1,lgtheta,lgtheta2,wave_angstrom,wave_cm,nu_Hz(iw),hm_bf_factor,stim_factor,&
                 ionization_factor,opacity,albedo) ! output: omega, and absorp_coeff
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
