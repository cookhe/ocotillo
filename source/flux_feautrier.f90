program flux_feautrier
!  
!                    _
!                 .'` '`.
!              .-", @ `, `.
!             '-=;      ;   `.
!                 \    :      `-.
!                 /    ';        `.
!                /      .'         `.
!                |     (      `.     `-.._
!                 \     \` ` `. \         `-.._
!                  `.   ;`-.._ `-`._.-. `-._   `-._
!                    `..'     `-.```.  `-._ `-.._.'
!                      `--..__..-`--'      `-.,'
!                         `._)`/
!                          /  /
!                         /--(
!                      -./,--'`-,
!                   ,^--(  
!                   ,--' `-,
!          ___ _ __   __ _ _ __ _ __ _____      __
!         / __| '_ \ / _` | '__| '__/ _ \ \ /\ / /
!         \__ \ |_) | (_| | |  | | | (_) \ V  V / 
!         |___/ .__/ \__,_|_|  |_|  \___/ \_/\_/  
!             | |                                 
!             |_|                                 
!
!
!  SPARROW - Spectral Post-processing AGN Radiaton transfeR On multiple Wavelengths
!
!  This program calculates 1D radiative transfer for simulation boxes, 
!  using the Feautrier method. The main application is for AGN disks; 
!  the temperatures and thus the opacities span several different regimes,
!  that are included in the code. 
!
!  Usage: setup the soft links, make, then run
!  
!             rt_setup 
!             make
!             ./source/flux_feautrier.x
!  
!  Authors: Harrison E. Cook
!           Wladimir Lyra  
!
!  © 2024
!
    use Auxiliary
    use Common
    use Disk
    use GasState
    use ContinuousOpacity
    use ReadAthena
    use FileIO
  
    implicit none

    type (pillar_case) :: p
!
! These quantities are the code output: mean intensity (U), flux (V) 
! and opacity (absorp_coeff).
!
    real, dimension(mz,nyloc,nxloc,nw) :: U
    real, dimension(nz,nyloc,nxloc,nw) :: V,absorp_coeff
!    
    real, dimension(nz,nyloc,nxloc) :: rho3d,temp3d
    real, dimension(mz) :: z
    real, dimension(nz) :: aa,bb,cc,dd
    real, dimension(nw) :: waves_angstrom
    real :: dz,z0,z1,dz1,dz2
    real :: start, finish
    real :: start_loop, finish_loop  
    real :: w0=3000,w1=5000
    real :: sigma_grey
    integer :: iw,ix,iy,iprocx,iprocy
    logical :: lgrey=.false.,lread_athena=.true.
    character(len=90) :: snapshot
    character(len=90) :: inputfile='./input.in'
!
    namelist /input/ z0,z1,w0,w1,sigma_grey,lgrey,lread_athena
!
! Start the time counter.
!
    call cpu_time(start)
!
! The code will take as argument a four string letter to mark the snapshots
! from Athena. If no argument is given, assume it is the 0000 snapshot.
! This is needed to launch the code on multiple snapshots from the command  
! line on a cluster for post-processing. 
!
    call getarg(1,snapshot)
    if (snapshot=='') snapshot="0000"
!  
!  Read the input namelist with the user-defined parameters.
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
! Read the user-defined values for temperature, density, and ionization threshold from input.in
!
   call read_temperature_input(inputfile)
   call read_density_input(inputfile)
   call read_gas_state_input(inputfile)
!
! Calculate the wavelength array. So far only allows linear.  
!   
   call calc_wavelength(w1,w0,waves_angstrom)
!
! Pre-calculate a number of wavelength-related but space- and time-independent quantities,
! to speed up the opacity calculation in run time. 
!
   call pre_calc_opacity_quantities(waves_angstrom)
! 
! The code primarily deals with athena output but also uses 
! artificial inputs for testing and benchmarking purposes if 
! lread_athena is false. 
!
   if (lread_athena) call read_athena_input(inputfile)
!
! Here starts the main loop through processors. The code will loop through each block of
! x and y processors, read the nxloc,nyloc,nzloc information, and go through the z processors
! to reconstitute a nxloc,nyloc,nz array. Then for cache-efficiency it will loop through 
! nxloc and nyloc and do operations on the z vertical 1D arrays (pillars).
!
   do iprocx=0,nprocx-1; do iprocy=0,nprocy-1   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  lroot = (iprocx==0) .and. (iprocy==0)
  if (lread_athena) then
    call read_from_athena(z,dz,rho3d,temp3d,iprocx,iprocy,snapshot)
  else
    call calc_grid(z1,z0,z,dz)
    call calc_density(p%rho,z)
    call calc_temperature(p%T,z)
  endif
  dz1=1./dz
  dz2=dz**2
!
!***********************************************************************
  !xy-dependent starts here
  if (lroot) call cpu_time(start_loop)
!***********************************************************************
  xloop: do ix=1,nxloc
    yloop: do iy=1,nyloc
      if (lread_athena) then
        p%rho =  rho3d(1:nz,iy,ix)
        p%T   = temp3d(1:nz,iy,ix)
      endif
      p%rho1=1./p%rho
      p%T1=1./p%T
      call calc_hydrogen_ion_frac(p)
      call solve_gas_state(p)
      p%electron_pressure = get_electron_pressure(p)
      p%e_scatter         = get_electron_thomson_scattering(p)
      p%theta             = 5040.* p%T1
      p%theta1            = 1./p%theta
      p%lgtheta           = log10(p%theta)
      p%lgtheta2          = p%lgtheta**2
      p%hm_bf_factor      = get_hydrogen_ion_bound_free(p)
!
! Loop over wavelengths
!
      wavelength: do iw=1,nw
        lfirst=lroot.and.(ix==1).and.(iy==1).and.(iw==1) 
!    
        if (lgrey) then
          call grey_parameters(p,sigma_grey)
        else
          p%source_function = get_source_function(p%T1,iw)
          p%stim_factor = get_hydrogen_stimulated_emission(iw,p%theta)
          call calc_opacity_and_albedo(p,iw)
        endif
!
! Populate coefficient arrays
!
        call fill_center_coeffs(aa,bb,cc,dd,p,dz2)
        call fill_boundary_coeffs(aa,bb,cc,dd,p,dz,dz2)
!
! Solve the system of equations
!
        if (lfirst) print*,'sum(a), sum(b), sum(c), sum(d)=',sum(aa), sum(bb), sum(cc), sum(dd)
        call tridag(aa, bb, cc, dd, U(n1:n2,iy,ix,iw))
        call update_ghosts(U(:,iy,ix,iw))
        call calc_flux(U(:,iy,ix,iw),p,dz1)
!
        absorp_coeff(:,iy,ix,iw)=p%opacity
        V(:,iy,ix,iw) = p%flux
!
        if (lfirst) then
           print*, 'min/max(U)', minval(U(n1:n2,iy,ix,iw)), maxval(U(n1:n2,iy,ix,iw))
           print*, 'min/max(V)', minval(V(:,iy,ix,iw)), maxval(V(:,iy,ix,iw))
        endif
!
      enddo wavelength
    enddo yloop
  enddo xloop
!
  if (lroot) call cpu_time(finish_loop)
!
! Calculate post-processing on one column for diagnostic: flux and intensities.
!
!  call calc_auxiliaries(U(:,ny,nx,:),absorp_coeff(:,ny,nx,:),dz,V,Ip,Im)
!
! Write output
!
  if (lroot) then
    call output_grid(z,waves_angstrom)
    !Output for diagnostic purposes 
    call output_ascii(U(:,nyloc,nxloc,:),absorp_coeff(:,nyloc,nxloc,:))
  endif
  call output_binary(U,absorp_coeff,V,iprocx,iprocy,snapshot)
  enddo;enddo
! 
! Finish the time counter and print the wall time
!
  call cpu_time(finish)
!
  print*,''
  print*,"Wall time = ",finish-start," seconds."
  print*,"Execution time =",(finish-start)/(1d0*nz*nyloc*nxloc*nw)*1e6,&
       " micro-seconds per frequency point per mesh point."
  print*,''
  print*,"Loop time = ",finish_loop-start_loop," seconds."
  print*,"Loop execution time =",(finish_loop-start_loop)*(nprocx*nprocy)/(1d0*nz*nyloc*nxloc*nw)*1e6,&
       " micro-seconds per frequency point per mesh point." 
  
endprogram flux_feautrier
!***********************************************************************
