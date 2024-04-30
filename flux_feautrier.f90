program flux_feautrier

  use auxiliary
  use grid

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

  integer :: greyindex,i,iwave,iz
  real :: sigma_sb=5.670374419d-5,sigma_grey=0.4
  real, dimension(nz) :: T,rho,B_grey,alpha_grey,kappa_grey,omega_grey

  real, dimension(nw) :: waves_cm,waves_angstrom
  real :: wave_cm,wave_angstrom

  real :: zeta,dz,z0,z1
  real :: rho0=1d-9,rho_floor=1d-24,H=1.496d14

  z0 = -4717816996570149.0
  z1 = 4717816996570149.0
  dz = (z1-z0)/nz
  do i=1,nz
   z(n1+i-1) = z0 + (i-1)*dz
   T(i) = 65000.
   rho(i) = rho0*exp(-.5*z(i-1+n1)**2/H**2)
   ! Implement a desity floor
   if (rho(i) < rho_floor) then 
      rho(i) = rho_floor
   endif ! density floor
  enddo ! populate z, T, rho arrays

!   rho=rho0*exp(-.5*z**2/H**2)
  overflow_limit = int(floor(log10(float_info_max)))
  
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
  do iwave=1,nw
     wave_cm = waves_cm(iwave)
     wave_angstrom = waves_angstrom(iwave)

     
     !damping_factor = h*c/wave_cm/kb/T
     !do i=1,nz
     !   if (damping_factor(i) > overflow_limit) then
     !      damping_factor(i) = overflow_limit
     !   endif
     !enddo
     !source_function(:,iwave) = 2*h*c**2/wave_cm**5 * 1/(exp(damping_factor)-1)
     
     !chi = 1.2398e4 / wave_angstrom
     !stim_factor = 1-10**(-chi*theta)
     
     !
     ! Do grey first 
     !
     
     absorp_coeff(:,iwave) = kappa_grey
     source_function(:,iwave) = B_grey
     omega(:,iwave) = omega_grey
                    
     ! Populate opacities at point i+1/2
     do iz=1, nz-1
       kappa_p(iz) = 0.5 * (absorp_coeff(iz+1,iwave) + absorp_coeff(iz,iwave))
     enddo
      
     ! Populate opacitiets at point i-1/2
     do iz=2, nz
       kappa_m(iz) = 0.5 * (absorp_coeff(iz,iwave) + absorp_coeff(iz-1,iwave))
     enddo

     ! Populate centers of arrays 
     do iz=2, nz-1
       aa(iz) = absorp_coeff(iz,iwave)**2 / kappa_m(iz)
       cc(iz) = absorp_coeff(iz,iwave)**2 / kappa_p(iz)
       zeta   = dz*dz * absorp_coeff(iz,iwave)**3 * (1 - omega(iz,iwave))
       bb(iz) = -(aa(iz) + cc(iz) + zeta)
       dd(iz) = -source_function(iz,iwave) * zeta
     enddo
   
     ! Populate boundary values
     aa(1) = 0
     zeta  = absorp_coeff(1,iwave) * dz**2 * (1 - omega(1,iwave)) / 4
     bb(1) = -(1/absorp_coeff(1,iwave) + dz + zeta) * absorp_coeff(1,iwave)**2
     cc(1) = 1/absorp_coeff(1,iwave) * absorp_coeff(1,iwave)**2
     dd(1) = -source_function(1,iwave) * zeta * absorp_coeff(1,iwave)**2
     
     aa(nz) = 1/absorp_coeff(nz,iwave) * absorp_coeff(nz,iwave)**2
     zeta   = absorp_coeff(nz,iwave) * dz**2 * (1 - omega(nz,iwave)) / 4
     bb(nz) = -(1/absorp_coeff(nz,iwave) + dz + zeta) * absorp_coeff(nz,iwave)**2
     cc(nz) = 0
     dd(nz) = -source_function(nz,iwave) * zeta * absorp_coeff(nz,iwave)**2
 
! solve the system of equations

     call tridag(aa, bb, cc, dd, U(n1:n2,iwave))
     call update_ghosts(U(:,iwave))
     
     call der(U(:,iwave),dU)
     V(:,iwave) = -1/absorp_coeff(:,iwave) * dU / dz

     Ip(:,iwave) = U(n1:n2,iwave) + V(:,iwave)
     Im(:,iwave) = U(n1:n2,iwave) - V(:,iwave)

     print*,'U(:,iwave)=',U(:,iwave)
     print*,''
     print*,'V(:,iwave)=',V(:,iwave)
     print*,''
     print*,'Ip(:,iwave)=',Ip(:,iwave)
     print*,''
     print*,'Im(:,iwave)=',Im(:,iwave)


  enddo ! wavelength loop
!
endprogram flux_feautrier
!********
