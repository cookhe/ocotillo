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
! These memory-intensive quantities are the code output: mean intensity (U), 
! flux (V) and opacity (absorp_coeff). Fortran is column major, but the 
! calculations will be done in z, so the data is assorted (z,y,x) to have 
! that direction be the fastest, for cache efficiency. 
!
    real, dimension(mz,nyloc,nxloc,nw) :: U
    real, dimension(nz,nyloc,nxloc,nw) :: V,absorp_coeff
!    
! Density and temperature are the athena input. The input will be by xy 
! processor and serial in the z direction.
!
    real, dimension(nz,nyloc,nxloc) :: rho3d,temp3d
    real, dimension(mz) :: z
!
    real, dimension(nz) :: a,b,c,d
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
!  Report which athena snapshot is being used.
!
    if (lread_athena) print*,"snapshot: ",snapshot
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
    procx: do iprocx=0,nprocx-1
      procy: do iprocy=0,nprocy-1   

        lroot = (iprocx==0) .and. (iprocy==0)
! 
! Read from athena input or construct a fake disk (density and temperature).
! The z direction is global, so only the root xy needs calculate it.
! All others can use it.
!
        if (lread_athena) then
          call read_from_athena(z,dz,rho3d,temp3d,iprocx,iprocy,snapshot)
        else
          if (lroot) call calc_grid(z1,z0,z,dz)
          call calc_density(p%rho,z)
          call calc_temperature(p%T,z)
        endif
!      
!  Shortcuts for the dz grid resolution.
!
        if (lroot) then 
          dz1=1./dz
          dz2=dz**2
        endif
!
! xy-dependent starts here.
!
        if (lroot) call cpu_time(start_loop)
!
        xloop: do ix=1,nxloc
          yloop: do iy=1,nyloc
!
! The fake disk has axisymmetric structure and thus the pillars for density
! and temperature in that case are defined already outside the xyloop. 
! Read here the athena input instead. 
!
            if (lread_athena) then
              p%rho =  rho3d(1:nz,iy,ix)
              p%T   = temp3d(1:nz,iy,ix)
            endif
!           
! Pre-define the pillars needed for the opacity calculation.
!            
            call wavelength_independent_pillars(p)
!
! Loop over wavelengths.
!
            wavelength: do iw=1,nw
              lfirst=lroot.and.(ix==1).and.(iy==1).and.(iw==1) 
!    
! These subroutines will calculate the opacity, albedo, and source function.  
!
              if (lgrey) then
                call grey_parameters(p,sigma_grey)
              else
                call calc_opacity_and_albedo(p,iw)
              endif
!
! Having the opacity, albedo, and source function, we can solve for the mean intensity. 
! This will be done via a tridiagonal system  
!
!    a_i U_{i-1} + b_i U_i + c_i U_{i+1} = d_i
!
! The next subroutine will calculate the coefficients a,b,c,d. 
!
              call get_tridag_coefficients(a,b,c,d,p,dz,dz2)
!
! Safety check for no NaNs in the coefficients. 
!
              if (lfirst) print*,'sum(a), sum(b), sum(c), sum(d)=',&
                                  sum(a), sum(b), sum(c), sum(d)
!
! Solve the tridiagonal system to get the mean intensity.
!
              call tridag(a,b,c,d,U(n1:n2,iy,ix,iw))
!              
! Calculate the flux V = -1./kappa  * dU/dz. 
! The mean intensity needs boundary conditions to calculate the derivative. 
!
              call update_ghosts(U(:,iy,ix,iw))
              call calc_flux(U(:,iy,ix,iw),p,dz1)
!
! Save opacity and flux into larger memomy arrays for output.
!
              absorp_coeff(:,iy,ix,iw)=p%opacity
              V(:,iy,ix,iw) = p%flux
!
! Safety check for no NaNs in mean intensity and flux. 
!
              if (lfirst) then
                print*, 'min/max(U)', minval(U(n1:n2,iy,ix,iw)), maxval(U(n1:n2,iy,ix,iw))
                print*, 'min/max(V)', minval(V(:,iy,ix,iw)), maxval(V(:,iy,ix,iw))
                print*, 'min/max(absorp_coeff)', minval(absorp_coeff(:,iy,ix,iw)), maxval(absorp_coeff(:,iy,ix,iw))
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
      enddo procy
    enddo procx
! 
! Finish the time counter and print the wall time
!
    call cpu_time(finish)
!
! Wall time only for the loop, in a single xy processor block. 
!
    print*,''
    print*,"processor column time = ",finish_loop-start_loop," seconds."
    print*,"Execution time =",(finish_loop-start_loop)/(1d0*nz*nyloc*nxloc*nw)*1e6,&
         " micro-seconds per wavelength point per mesh point."
!
! Wall time for the whole code
!
    print*,''
    print*,"Wall time = ",finish-start," seconds."
    print*,"Execution time =",(finish-start)/(1d0*nz*ny*nx*nw)*1e6,&
         " micro-seconds per wavelength point per mesh point."
!
endprogram flux_feautrier
!***********************************************************************
