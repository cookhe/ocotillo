program flux_feautrier

  use auxiliary
  use grid
  use disk, only:temperature_gaussian

  implicit none

  real, dimension(mz) :: z

  real, dimension(mz,nw) :: U
  real, dimension(nz,nw) :: V,Ip,Im
  !real, dimension(nz,nw) :: kappa_H_bf,kappa_H_ff,kappa_Hm_bf,kappa_Hm_ff,kappa_rad
  real, dimension(nz,nw) :: absorp_coeff,source_function,omega

  real, dimension(nz) :: aa,bb,cc,dd,kappa_p,kappa_m,dU
  
  !real :: alpha_e = 0.6648e-24 ! coefficient
  integer :: overflow_limit
  real :: float_info_max=3.9085d307

  integer :: greyindex,i,iw,iz
  real :: sigma_sb=5.670374419d-5,sigma_grey=0.4
  real, dimension(nz) :: T,rho,B_grey,alpha_grey,kappa_grey,omega_grey

  real, dimension(nw) :: waves_cm,waves_angstrom
  real :: wave_cm,wave_angstrom

  real :: zeta,dz,z0,z1,rho0,rho_floor,H

  real :: start, finish
!
  namelist /input/ z0,z1,rho0,rho_floor,H
!  
!  Read the input namelist with the user-defined parameters. 
!  
  open(20,file='input.in')
  read(20,nml=input)
  close(20)
!    
  dz = (z1-z0)/nz
  arrays: do i=1,nz
   z(n1+i-1) = z0 + (i-1)*dz
   T(i) = 65000.
   rho(i) = rho0*exp(-.5*z(i-1+n1)**2/H**2)
   ! Implement a desity floor
   density_floor: if (rho(i) < rho_floor) then 
      rho(i) = rho_floor
   endif density_floor
  enddo arrays

  overflow_limit = int(floor(log10(float_info_max)))
  
  ! Apply temperature profile
  call temperature_gaussian(T,z)
  print*,"maxval(T), minval(T)",maxval(T), minval(T)

  !n = rho*mp1
  !call hydrogen_ionization_fraction(rho,T,NHII_NHINHII)

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
  
  greyindex = 0
  B_grey = sigma_sb*T**4
  alpha_grey = 3.68e22 * T**(-3.5) *  rho
  kappa_grey = (alpha_grey+sigma_grey)*rho
  omega_grey = sigma_grey / (alpha_grey + sigma_grey)

  ! Loop over wavelength

  call cpu_time(start)

  wavelength: do iw=1,nw
     wave_cm = waves_cm(iw)
     wave_angstrom = waves_angstrom(iw)

     
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
     ! Do grey first 
     !
     
     absorp_coeff(:,iw) = kappa_grey
     source_function(:,iw) = B_grey
     omega(:,iw) = omega_grey
                    
     ! Populate opacities at point i+1/2
     do iz=1, nz-1
       kappa_p(iz) = 0.5 * (absorp_coeff(iz+1,iw) + absorp_coeff(iz,iw))
     enddo
      
     ! Populate opacitiets at point i-1/2
     do iz=2, nz
       kappa_m(iz) = 0.5 * (absorp_coeff(iz,iw) + absorp_coeff(iz-1,iw))
     enddo

     ! Populate centers of arrays 
     do iz=2, nz-1
       aa(iz) = absorp_coeff(iz,iw)**2 / kappa_m(iz)
       cc(iz) = absorp_coeff(iz,iw)**2 / kappa_p(iz)
       zeta   = dz*dz * absorp_coeff(iz,iw)**3 * (1 - omega(iz,iw))
       bb(iz) = -(aa(iz) + cc(iz) + zeta)
       dd(iz) = -source_function(iz,iw) * zeta
     enddo
   
     ! Populate boundary values
     aa(1) = 0.
     zeta  = absorp_coeff(1,iw) * dz**2 * (1 - omega(1,iw)) / 4
     bb(1) = -(1/absorp_coeff(1,iw) + dz + zeta) * absorp_coeff(1,iw)**2
     cc(1) = 1/absorp_coeff(1,iw) * absorp_coeff(1,iw)**2
     dd(1) = -source_function(1,iw) * zeta * absorp_coeff(1,iw)**2
     
     aa(nz) = 1/absorp_coeff(nz,iw) * absorp_coeff(nz,iw)**2
     zeta   = absorp_coeff(nz,iw) * dz**2 * (1 - omega(nz,iw)) / 4
     bb(nz) = -(1/absorp_coeff(nz,iw) + dz + zeta) * absorp_coeff(nz,iw)**2
     cc(nz) = 0.
     dd(nz) = -source_function(nz,iw) * zeta * absorp_coeff(nz,iw)**2
 
! solve the system of equations

     call tridag(aa, bb, cc, dd, U(n1:n2,iw))

  enddo wavelength
!
  call cpu_time(finish)
  print*,"Wall time = ",finish-start," seconds."

  !put this inside a subroutine, like, calculate_auxiliaries

  do iw=1,nw
    call update_ghosts(U(:,iw))
    call der(U(:,iw),dU)
    V(:,iw) = -1/absorp_coeff(:,iw) * dU / dz

    Ip(:,iw) = U(n1:n2,iw) + V(:,iw)
    Im(:,iw) = U(n1:n2,iw) - V(:,iw)
  enddo
!
  open(10,file="intensity.dat",status="replace",action='write')
  do i=1,nz
    do iw=1,nw
      write(unit=10,FMT=*) i,iw,U(n1+i-1,iw),V(i,iw),Ip(i,iw),Im(i,iw)
    enddo
  enddo
  close(10)
  
endprogram flux_feautrier
!********
