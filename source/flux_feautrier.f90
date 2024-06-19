program flux_feautrier

  use Auxiliary
  use Common
  use Disk
  use GasState
  use ContinuousOpacity
  use Grey
  use ReadC

  implicit none

  real, dimension(mz,nw) :: U
  real, dimension(nz,nw) :: V,Ip,Im
  real, dimension(nz,nw) :: absorp_coeff
  real, dimension(mz) :: z
  real, dimension(nz) :: aa,bb,cc,dd

  integer :: iw,ix,iy
  integer :: log_overflow_limit

  real :: sigma_grey
  real, dimension(nz,ny,nx) :: rho3d
  real, dimension(nz) :: rho,T
  real, dimension(nz) :: opacity, albedo
  real, dimension(nw) :: waves_cm,waves_angstrom

  real, dimension(nz) :: NHII_NHINHII,n,nHII,nHI,ne,ionization_factor
  real, dimension(nz) :: e_scatter,theta,electron_pressure,hm_bf_factor
  real, dimension(nz) :: stim_factor,source_function

  real :: wave_cm,wave_angstrom
  real :: dz,z0,z1
  real :: start, finish
  real :: w0=3000,w1=5000
!
  logical :: lgrey=.false.,lread_athena=.true.
!
  !real, dimension(nz,nw) :: kappa_H_bf,kappa_H_ff,kappa_Hm_bf,kappa_Hm_ff,kappa_rad
  !real, dimension(nz) :: kappa_p,kappa_m 
!
  namelist /input/ z0,z1,w0,w1,sigma_grey,lgrey,lread_athena
!  
!  Read the input namelist with the user-defined parameters. 
!  
  open(20,file='./input.in')
  read(20,nml=input)
  close(20)
!    
  log_overflow_limit = int(floor(log10(float_info_max)))
!
! Calculate the grid variables
!
  call calc_wavelength(w1,w0,waves_angstrom,waves_cm)
  if (.not.lread_athena) then 
     call calc_grid(z1,z0,z,dz)
     call calc_density(rho,z)
  else
     call read_from_athena(z,dz,rho3d)
  endif

  call calc_temperature(T,z)
!
! Start the counter
!
      call cpu_time(start)    
!*******************************
!xy-dependent starts here
!*******************************
  xloop: do ix=1,nx
    yloop: do iy=1,ny
      if (lread_athena) rho = rho3d(1:nz,iy,ix)
      call calc_hydrogen_ion_frac(rho,T,NHII_NHINHII)
      call solve_gas_state(rho,NHII_NHINHII,n,nHI,nHII,ne,ionization_factor)
      call calc_electron_pressure(ne,T,electron_pressure)
      e_scatter    = get_electron_thomson_scattering(n,nHII)
      theta        = get_theta(T)
      hm_bf_factor = get_hydrogen_ion_bound_free(electron_pressure,theta)
!
! Do a sanity check for grey RT
!
      if (lgrey.and.(nw /= 1)) then
         print*,"For Grey RT use only one wavelength. Switch nw=1 in commons.f90"
         stop
      endif
!
! Loop over wavelengths
!
      wavelength: do iw=1,nw
         lfirst=(ix==1).and.(iy==1).and.(iw==1) 
         wave_cm = waves_cm(iw)
         if (lfirst) print*, 'waves_cm',waves_cm
         wave_angstrom = waves_angstrom(iw)
!    
         if (lgrey) then
            call grey_parameters(rho,T,sigma_grey,source_function,opacity,albedo)
         else
            call calc_source_function(wave_cm,T,source_function,log_overflow_limit) !output: source_function
            call calc_hydrogen_stimulated_emission(wave_angstrom,theta,stim_factor) !output: stim_factor
            call calc_opacity_and_albedo(e_scatter,rho,ne,NHII_NHINHII,nHI,nHII,&
                 T,wave_angstrom,hm_bf_factor,stim_factor,ionization_factor,opacity,albedo) ! output: omega, and absorp_coeff
         endif
!
! Populate coefficient arrays
!
         call fill_center_coeffs(aa,bb,cc,dd,opacity,albedo,source_function,dz)
         call fill_boundary_coeffs(aa,bb,cc,dd,opacity,albedo,source_function,dz)    
!
! Solve the system of equations
!
         call tridag(aa, bb, cc, dd, U(n1:n2,iw))
         if (lfirst) print*, 'min/max(U)', minval(U(n1:n2,iw)), maxval(U(n1:n2,iw))
!
         absorp_coeff(:,iw)=opacity
!
      enddo wavelength
    enddo yloop
  enddo xloop
!
! Finish the counter and print the wall time
!
  call cpu_time(finish)
  print*,"Wall time = ",finish-start," seconds."
  print*,"Execution time =",(finish-start)/(nz*ny*nx*nw)*1e6,&
       " micro-seconds per frequency point per mesh point." 
!
! Calculate post-processing quantities: flux and intensities. 
!
  call calc_auxiliaries(U,absorp_coeff,dz,V,Ip,Im)
!
! Write output
!
  call output_grid(z)
  call output_data(U,V,Ip,Im)
! 
endprogram flux_feautrier
!********
