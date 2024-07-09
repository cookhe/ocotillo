module Auxiliary

  use Common

  implicit none
  private
  
  public :: tridag,update_ghosts,der
  public :: calc_auxiliaries,calc_flux
  public :: fill_center_coeffs,fill_boundary_coeffs

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
  subroutine calc_flux(U,p,dz1)

    real, dimension(mz) :: U
    real, dimension(nz) :: dU
    type (pillar_case) :: p
    real :: dz1

    call der(U,dU)
    p%flux = -1*p%opacity1*dU*dz1

  endsubroutine calc_flux
!************************************************************************************
  subroutine fill_center_coeffs(aa,bb,cc,dd,p,dz2)
    real, dimension(nz), intent(inout) :: aa,bb,cc,dd
    real, dimension(nz) :: kappa_m,kappa_p
    real, intent(in) :: dz2
    real :: zeta
    integer :: iz
    type (pillar_case) :: p
  
    ! Populate opacities at point i+1/2
    do iz=1, nz-1
      kappa_p(iz) = 0.5 * (p%opacity(iz+1) + p%opacity(iz))
    enddo
     
    ! Populate opacities at point i-1/2
    do iz=2, nz
      kappa_m(iz) = 0.5 * (p%opacity(iz) + p%opacity(iz-1))
    enddo

    ! Populate centers of arrays 
    do iz=2, nz-1
      aa(iz) = p%opacity2(iz) / kappa_m(iz)
      cc(iz) = p%opacity2(iz) / kappa_p(iz)
      zeta   = dz2 * p%opacity(iz)**3 * (1 - p%albedo(iz))
      bb(iz) = -(aa(iz) + cc(iz) + zeta)
      dd(iz) = -p%source_function(iz) * zeta
    enddo

    if (lfirst) print*, ' Filled center coeffs'
    
  endsubroutine fill_center_coeffs
!************************************************************************************
  subroutine fill_boundary_coeffs(aa,bb,cc,dd,p,dz,dz2)
    
    real, dimension(nz), intent(inout) :: aa,bb,cc,dd
    real, intent(in) :: dz,dz2
    real :: zeta

    type (pillar_case) :: p
  
    ! Populate boundary values
    aa(1) = 0.
    zeta  = p%opacity(1) * dz2 * (1 - p%albedo(1))*.25
    bb(1) = -(p%opacity1(1) + dz + zeta) * p%opacity2(1)
    cc(1) = p%opacity(1)
    dd(1) = -p%source_function(1) * zeta * p%opacity2(1)
    
    aa(nz) = p%opacity(nz)
    zeta   = p%opacity(nz) * dz2 * (1 - p%albedo(nz))*.25
    bb(nz) = -(p%opacity1(nz) + dz + zeta) * p%opacity2(nz)
    cc(nz) = 0.
    dd(nz) = -p%source_function(nz) * zeta * p%opacity2(nz)

    if (lfirst) print*, ' Filled boundary coeffs'
    
  endsubroutine fill_boundary_coeffs
!************************************************************************************
endmodule Auxiliary
