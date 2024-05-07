program flux_feautrier

  use auxiliary
  use common
  use disk
  use gas_state

  implicit none

  real, dimension(mz,nw) :: U
  real, dimension(nz,nw) :: V,Ip,Im
  real, dimension(nz,nw) :: absorp_coeff,source_function,omega
  real, dimension(mz) :: z
  real, dimension(nz) :: aa,bb,cc,dd

  integer :: iw
  integer :: overflow_limit

  real :: sigma_grey
  real, dimension(nz) :: T,rho,n,B_grey,alpha_grey,kappa_grey,omega_grey
  real, dimension(nw) :: waves_cm,waves_angstrom

  real, dimension(nz) :: NHII_NHINHII

  real :: wave_cm,wave_angstrom
  real :: dz,z0,z1
  real :: start, finish
!
  !real, dimension(nz,nw) :: kappa_H_bf,kappa_H_ff,kappa_Hm_bf,kappa_Hm_ff,kappa_rad
  !real, dimension(nz) :: kappa_p,kappa_m 
  !real :: alpha_e = 0.6648e-24 ! coefficient
!
  namelist /input/ z0,z1,sigma_grey
!  
!  Read the input namelist with the user-defined parameters. 
!  
  open(20,file='input.in')
  read(20,nml=input)
  close(20)
!    
  overflow_limit = int(floor(log10(float_info_max)))
!
! Calculate the grid variables
!
  call calc_grid(z1,z0,z,dz)
  call calc_density(rho,z)
  ! print*, ' rho=', rho
  ! rho = merge(1d-10,rho,rho<1d-10)
  ! print*, ' rho=', rho
  call calc_temperature(T,z)
  call hydrogen_ion_frac(rho,T,NHII_NHINHII)
  print*,NHII_NHINHII
  
  ! n = rho*mp1
  !nHII = NHII_NHINHII * n
  !nHI = n - nHII
  !ne = nHII

  !do i=1,nz
  !  if (nHI(i) /= 0) then
  !    ionization_factor[i] = 1/(1 + nHII[i]/nHI[i])
  !  endif
  !enddo
  
  !e_scatter = alpha_e * (nHII / n)  
  !theta = 5040./T
  !electron_pressure = ne * kb * T
  !hm_bf_factor = 4.158e-10 * electron_pressure * theta**(5./2) * 10**(0.754 * theta)
  
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
    wave_angstrom = waves_angstrom(iw)
!    
    !damping_factor = h*c/wave_cm/kb/T
    !do i=1,nz
    !   if (damping_factor(i) > overflow_limit) then
    !      damping_factor(i) = overflow_limit
    !   endif
    !enddo
    !source_function(:,iw) = 2*h*c**2/wave_cm**5 * 1/(exp(damping_factor)-1)
    
    !chi = 1.2398e4 / wave_angstrom
    !stim_factor = 1-10**(-chi*theta)
!
! Populate coefficient arrays
!
    call fill_center_coeffs(aa,bb,cc,dd,absorp_coeff(:,iw),omega(:,iw),source_function(:,iw),dz)
    call fill_boundary_coeffs(aa,bb,cc,dd,absorp_coeff(:,iw),omega(:,iw),source_function(:,iw),dz)    
!
! Solve the system of equations
!
    call tridag(aa, bb, cc, dd, U(n1:n2,iw))
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
