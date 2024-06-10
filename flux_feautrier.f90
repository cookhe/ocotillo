program flux_feautrier

  use Auxiliary
  use Common
  use Disk
  use GasState

  implicit none

  real, dimension(mz,nw) :: U
  real, dimension(nz,nw) :: V,Ip,Im
  real, dimension(nz,nw) :: absorp_coeff,source_function,omega
  real, dimension(mz) :: z
  real, dimension(nz) :: aa,bb,cc,dd

  integer :: iw,i
  integer :: log_overflow_limit

  real :: sigma_grey
  real, dimension(nz) :: T,rho,B_grey,alpha_grey,kappa_grey,omega_grey
  real, dimension(nw) :: waves_cm,waves_angstrom

  real, dimension(nz) :: NHII_NHINHII,n,nHII,nHI,ne,ionization_factor
  real, dimension(nz) :: e_scatter,theta,electron_pressure,hm_bf_factor

  real :: wave_cm,wave_angstrom
  real :: dz,z0,z1
  real :: start, finish
!
  !real, dimension(nz,nw) :: kappa_H_bf,kappa_H_ff,kappa_Hm_bf,kappa_Hm_ff,kappa_rad
  !real, dimension(nz) :: kappa_p,kappa_m 
  real :: alpha_e = 0.6648e-24 ! coefficient
!
  namelist /input/ z0,z1,sigma_grey
!  
!  Read the input namelist with the user-defined parameters. 
!  
  open(20,file='input.in')
  read(20,nml=input)
  close(20)
!    
  log_overflow_limit = int(floor(log10(float_info_max)))
!
! Calculate the grid variables
!
  do iw=1,nw
    waves_angstrom(iw) = w1 + ((iw-1)*dw)
  enddo
  print*, 'waves_angstrom', waves_angstrom
  waves_cm = waves_angstrom*1d-8
  call calc_grid(z1,z0,z,dz)
  call calc_density(rho,z)
  call calc_temperature(T,z)
  call hydrogen_ion_frac(rho,T,NHII_NHINHII)
  call solve_gas_state(rho,NHII_NHINHII,n,nHI,nHII,ne,ionization_factor)
  
  e_scatter = alpha_e * (nHII / n)  
  theta = 5040./T
  electron_pressure = ne * k_cgs * T
  hm_bf_factor = 4.158e-10 * electron_pressure * theta**(5./2) * 10**(0.754 * theta)
  
  B_grey = sigma_sb*T**4
  alpha_grey = 3.68e22 * T**(-3.5) *  rho
  kappa_grey = (alpha_grey+sigma_grey)*rho
  omega_grey = sigma_grey / (alpha_grey + sigma_grey)
!
! Start the counter
!
  call cpu_time(start)  
!
! Do grey RT first 
!
  absorp_coeff(:,igrey) = kappa_grey
  source_function(:,igrey) = B_grey
  omega(:,igrey) = omega_grey
!
! Loop over wavelengths
!
  wavelength: do iw=1,nw
    wave_cm = waves_cm(iw)
    ! print*, 'waves_cm',waves_cm
    wave_angstrom = waves_angstrom(iw)
!    
    call calc_source_function(wave_cm,T,source_function,iw,log_overflow_limit) !output: source_function
!     
    !call calc_albedo_and_opacity() ! output: omega, and absorp_coeff

    ! next call kappa_rad, which will be used for albedo and absorption
    
    !chi = 1.2398e4 / wave_angstrom
    !stim_factor = 1-10**(-chi*theta)
    !
    !kappa_H_bf[:,iwave]  = specb.xsec_bfHI(np.array([wave_angstrom]), T, m=6)
    !kappa_H_bf[:,iwave] *= stim_factor * ionization_factor
    !kappa_H_ff[:,iwave]   = specb.xsec_ffH(np.array(wave_angstrom), T)
    !kappa_H_ff[:,iwave]  *= stim_factor * ionization_factor
    !
    !alpha_Hmbf          = specb.xsec_bfHminus(np.array([wave_angstrom]))
    !kappa_Hm_bf[:,iwave]  = hm_bf_factor * alpha_Hmbf
    !kappa_Hm_bf[:,iwave] *= stim_factor * ionization_factor
    !
    !kappa_Hm_ff[:,iwave]  = specb.xsec_ffHminus(ne, np.array([wave_angstrom]), T)
    !kappa_Hm_ff[:,iwave] *= ionization_factor

    !kappa_rad[:,iwave] = kappa_H_bf[:,iwave] + kappa_H_ff[:,iwave] + kappa_Hm_bf[:,iwave] + kappa_Hm_ff[:,iwave]
    !kappa_rad[:,iwave] += e_scatter

    ! placeholder until albedo properly calculated
    omega(:,iw) = omega_grey
    !omega[:,iwave] = (e_scatter) / (kappa_rad[:,iwave])
    !kappa_rad[:,iwave] /= mp
    
    ! placeholder until absorption coefficient properly calculated
    absorp_coeff(:,iw) = kappa_grey
    !absorp_coeff[:,iwave] = kappa_rad[:,iwave] * rho
       
!
! Populate coefficient arrays
!
    call fill_center_coeffs(aa,bb,cc,dd,absorp_coeff(:,iw),omega(:,iw),source_function(:,iw),dz)
    call fill_boundary_coeffs(aa,bb,cc,dd,absorp_coeff(:,iw),omega(:,iw),source_function(:,iw),dz)    
!
! Solve the system of equations
!
    call tridag(aa, bb, cc, dd, U(n1:n2,iw))
    print*, 'min/max(U)', minval(U(n1:n2,iw)), maxval(U(n1:n2,iw))
    ! print*, source_function(:,iw)
!
  enddo wavelength
!
! Finish the counter and print the wall time
!
  call cpu_time(finish)
  print*,"Wall time = ",finish-start," seconds."
!
! Calculate post-processing quantities: flux and intensities. 
!
  call calc_auxiliaries(U,absorp_coeff,dz,V,Ip,Im)
!
! Write output
!
  call output_data(U,V,Ip,Im)
! 
endprogram flux_feautrier
!********
