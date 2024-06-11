module Auxiliary

  use Common

  implicit none
  private
  
  public :: tridag,update_ghosts,der
  public :: calc_auxiliaries
  public :: fill_center_coeffs,fill_boundary_coeffs
  public :: output_data
  public :: calc_source_function

contains
!************************************************************************************
  subroutine tridag(a,b,c,r,u)
!
!  Solves a tridiagonal system.
!  Imported from numerical recipes.    
!
    real, dimension(:), intent(in) :: a,b,c,r
    real, dimension(:), intent(out) :: u
    real, dimension(size(b)) :: gam
    integer :: n,j
    real :: bet
!
    n=size(b)
    bet=b(1)
    if (bet==0.0) then
      print*,"ERROR tridiag stage 1: bet=b(1)) = ", bet
      stop
    endif
!
    u(1)=r(1)/bet
    do j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      if (bet==0.0) then
        print*,"ERROR tridiag stage 2: bet = b(j)-a(j)*gam(j) =", bet
        stop
      endif
      u(j)=(r(j)-a(j)*u(j-1))/bet
    enddo
!
    do j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
    enddo
!
  endsubroutine tridag
!************************************************************************************
  subroutine update_ghosts(f)
!
!  Update the ghost zones of an array
!  using constant gradient. Used on U.     
!
    real, dimension(:), intent(inout) :: f
    integer :: i
!      
    do i=1,ng
      f(n1-i)=2*f(n1) - f(n1+i)
      f(n2+i)=2*f(n2) - f(n2-i)
    enddo
!
  endsubroutine update_ghosts
!************************************************************************************
  subroutine der(f,df)
!
!  Sixth-order first derivative.
!    
    real, dimension(mz) :: f
    real, dimension(nz) :: df
!      
    intent(in) :: f
    intent(out) :: df
!
    df=1./60.*(+ 45.0*(f(n1+1:n2+1)-f(n1-1:n2-1)) &
               -  9.0*(f(n1+2:n2+2)-f(n1-2:n2-2)) &
               +      (f(n1+3:n2+3)-f(n1-3:n2-3)))
!    
  endsubroutine der
!************************************************************************************
  subroutine calc_auxiliaries(U,absorp_coeff,dz,V,Ip,Im)

    real, dimension(mz,nw) :: U
    real, dimension(nz,nw) :: V,Ip,Im,absorp_coeff
    real, dimension(nz) :: dU
    real :: dz
    integer :: iw
!
    do iw=1,nw
      call update_ghosts(U(:,iw))
      call der(U(:,iw),dU)
      V(:,iw) = -1/absorp_coeff(:,iw) * dU / dz
!
      Ip(:,iw) = U(n1:n2,iw) + V(:,iw)
      Im(:,iw) = U(n1:n2,iw) - V(:,iw)
    enddo
!
  endsubroutine calc_auxiliaries
!************************************************************************************
  subroutine fill_center_coeffs(aa,bb,cc,dd,absorp_coeff,omega,source_function,dz)
    real, dimension(nz), intent(inout) :: aa,bb,cc,dd
    real, dimension(nz), intent(in) :: absorp_coeff,omega,source_function
    real, dimension(nz) :: kappa_m,kappa_p
    real, intent(in) :: dz
    real :: zeta
    integer :: iz
  
    ! Populate opacities at point i+1/2
    do iz=1, nz-1
      kappa_p(iz) = 0.5 * (absorp_coeff(iz+1) + absorp_coeff(iz))
    enddo
     
    ! Populate opacities at point i-1/2
    do iz=2, nz
      kappa_m(iz) = 0.5 * (absorp_coeff(iz) + absorp_coeff(iz-1))
    enddo

    ! Populate centers of arrays 
    do iz=2, nz-1
      aa(iz) = absorp_coeff(iz)**2 / kappa_m(iz)
      cc(iz) = absorp_coeff(iz)**2 / kappa_p(iz)
      zeta   = dz*dz * absorp_coeff(iz)**3 * (1 - omega(iz))
      bb(iz) = -(aa(iz) + cc(iz) + zeta)
      dd(iz) = -source_function(iz) * zeta
    enddo

    print*, ' Filled center coeffs'
    
  endsubroutine fill_center_coeffs
!************************************************************************************
  subroutine fill_boundary_coeffs(aa,bb,cc,dd,absorp_coeff,omega,source_function,dz)
    
    real, dimension(nz), intent(inout) :: aa,bb,cc,dd
    real, dimension(nz), intent(in) :: absorp_coeff,omega,source_function
    real, intent(in) :: dz
    real :: zeta
  
    ! Populate boundary values
    aa(1) = 0.
    zeta  = absorp_coeff(1) * dz**2 * (1 - omega(1)) / 4
    bb(1) = -(1/absorp_coeff(1) + dz + zeta) * absorp_coeff(1)**2
    cc(1) = 1/absorp_coeff(1) * absorp_coeff(1)**2
    dd(1) = -source_function(1) * zeta * absorp_coeff(1)**2
    
    aa(nz) = 1/absorp_coeff(nz) * absorp_coeff(nz)**2
    zeta   = absorp_coeff(nz) * dz**2 * (1 - omega(nz)) / 4
    bb(nz) = -(1/absorp_coeff(nz) + dz + zeta) * absorp_coeff(nz)**2
    cc(nz) = 0.
    dd(nz) = -source_function(nz) * zeta * absorp_coeff(nz)**2

    print*, ' Filled boundary coeffs'
    
  endsubroutine fill_boundary_coeffs
!************************************************************************************
  subroutine output_data(U,V,Ip,Im)

    real, dimension(mz,nw) :: U
    real, dimension(nz,nw) :: V,Ip,Im
    integer :: i,iw
!
    open(10,file="intensity.dat",status="replace",action='write')
    do i=1,nz
      do iw=1,nw
        write(unit=10,FMT=*) i,iw,U(n1+i-1,iw),V(i,iw),Ip(i,iw),Im(i,iw)
      enddo
    enddo
    close(10)

  endsubroutine output_data
!************************************************************************************
  subroutine calc_source_function(wave_cm,T,source_function,log_overflow_limit)

    integer :: i,log_overflow_limit
    real, intent(in) :: wave_cm
    real, dimension(nz) :: T,damping_factor
    real, dimension(nz) :: source_function
    
    damping_factor = h_planck*c_light_cgs/(wave_cm*k_cgs*T)
    do i=1,nz
      if (damping_factor(i) > log_overflow_limit) then
        damping_factor(i) = log_overflow_limit
      endif
    enddo
    
    source_function = 2*h_planck*c_light_cgs**2/wave_cm**5 * 1/(exp(damping_factor)-1)
    
  endsubroutine calc_source_function
!************************************************************************************    
  endmodule Auxiliary
