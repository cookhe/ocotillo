program flux_feautrier

  integer, parameter :: nw=1
  integer, parameter :: nz=640
  integer, parameter :: ng=3

  real, dimension(nz) :: z

  real, dimension(nz+2*ng,nw) :: U
  real, dimension(nz,nw) :: V,Ip,Im
  !real, dimension(nz,nw) :: kappa_H_bf,kappa_H_ff,kappa_Hm_bf,kappa_Hm_ff,kappa_rad
  real, dimension(nz,nw) :: absorp_coeff,source_function,omega

  real, dimension(nz) :: aa,bb,cc,dd
  
  !real :: alpha_e = 0.6648e-24 ! coefficient
  integer :: overflow_limit_exp
  !real :: float_info_max = 1.797693134862316E+308

  integer :: greyindex
  real :: sigma_sb=5.670374419d-5,sigma_grey=0.4
  real, dimension(nz) :: T,rho,B_grey,alpha_grey,kappa_grey,omega_grey,dU

  real, dimension(nw) :: waves_cm,waves_angstrom
  real :: wave_cm,wave_angstrom

  real :: zeta,dz,z0,z1
  real :: H=1.496d14,rho0=1d-9

  z0 = -4717816996570149.0
  z1 = 4717816996570149.0
  dz = (z1-z0)/nz
  do i=1,nz
     z[i] = z0 + (i-1)*dz
     T[i] = 65000.
  enddo
  
  rho=rho0*exp(-.5*z**2/H**2)
  overflow_limit = int(floor(log10(float_info)))
  
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
     !call TDMAsolver(aa, bb, cc, dd, U(:,iwave))
     !call update_ghosts(U(:,iwave))
     !call der6(U(:,iwave),dU)
     V(:,iwave) = -1/absorp_coeff(:,iwave) * dU / dz

     Ip(:,iwave) = U(l1:l2,iwave) + V(:,iwave)
     Im(:,iwave) = U(l1:l2,iwave) - V(:,iwave)

  enddo


!
!subroutine der6(f,df)
!!
!  real, dimension (mx), intent(in)  :: f
!  real, dimension (nx), intent(out) :: df
!!                                                                                                                                                                                !  
!  real, dimension (nx) :: fac
!!                                                                                                                                                                                ! 
!  fac=(1./60)*dx_1(l1:l2)
!  df=fac*(+ 45.0*(f(l1+1:l2+1)-f(l1-1:l2-1)) &
!       -  9.0*(f(l1+2:l2+2)-f(l1-2:l2-2)) &
!       +      (f(l1+3:l2+3)-f(l1-3:l2-3)))
!!                                                                                                                                                                                !  
!    endsubroutine der_x

endprogram flux_feautrier
