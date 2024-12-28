module Auxiliary
  !! Auxiliary subroutines used throughout the code.

  use Common

  implicit none
  private
  
  public :: tridag,update_ghosts,der
  public :: calc_intensity,calc_flux
  public :: get_tridag_coefficients

contains
!====================================================================================
  subroutine tridag(a,b,c,r,u)
    !! Solves a tridiagonal system of equations.
    !! Imported from numerical recipes.
!
    real, dimension(:), intent(in) :: a,b,c,r
    !! Coefficients for the left-hand side (a, b, c) 
    !! and value of the right-hand side (r)
    real, dimension(:), intent(out) :: u
    !! Value for which to solve.
    real, dimension(size(b)) :: gam
    !! Intermediate step variable.
    integer :: n,j
    !! n : size of b, j : do-loop integer
    real :: bet
    !! Intermediate step variable.
!
    n=size(b)
    bet=b(1)
    if (bet==0.0) then
      print*,"ERROR tridiag stage 1: bet=b(1)) = ", bet, ": cannot be zero."
      stop
    endif
!
    u(1)=r(1)/bet
    do j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      if (bet==0.0) then
        print*,"ERROR tridiag stage 2: bet = b(j)-a(j)*gam(j) =", bet, ": cannot be zero."
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
!====================================================================================
  subroutine update_ghosts(f)
    !! Update the ghost zones of an array 
    !! using constant gradient. Used on U.     
!
    real, dimension(:), intent(inout) :: f
    !! Array to update.
    integer :: i
    !! Looping integer.
!      
    do i=1,ng
      f(n1-i)=2*f(n1) - f(n1+i)
      f(n2+i)=2*f(n2) - f(n2-i)
    enddo
!
  endsubroutine update_ghosts
!====================================================================================
  subroutine der(f,df)
    !! Calculate the sixth-order first derivative of an array.
!    
    real, dimension(mz) :: f
    !! Array to calculate the derivative including ghost zones.
    real, dimension(nz) :: df
    !! Derivative of the array.

!      
    intent(in) :: f
    intent(out) :: df
!
    df=1./60.*(+ 45.0*(f(n1+1:n2+1)-f(n1-1:n2-1)) &
               -  9.0*(f(n1+2:n2+2)-f(n1-2:n2-2)) &
               +      (f(n1+3:n2+3)-f(n1-3:n2-3)))
!    
  endsubroutine der
!====================================================================================
  subroutine calc_intensity(U,V,Ip,Im)
    !* Calculate upwards (Ip) and downwards (Im) intensities
    !  from U and V.
    ! 
    !  $$ I_{\rm p} = U + V $$
    !* $$ I_{\rm m} = U - V $$
!
    real, dimension(mz,nw) :: U
    !! Mean intensity array.
    real, dimension(nz,nw) :: V,Ip,Im
    !! Flux, upward intensity, and downward intensity.
    integer :: iw
    !! Wavelength loop index.
!
    do iw=1,nw
      Ip(:,iw) = U(n1:n2,iw) + V(:,iw)
      Im(:,iw) = U(n1:n2,iw) - V(:,iw)
    enddo
!
  endsubroutine calc_intensity
!====================================================================================
  subroutine calc_flux(U,p,dz1)
    !! Calculate flux (V) from the mean intensity (U) and opacity (\(\kappa\))
    !! $$ V = -\frac{1}{\kappa} * \frac{dU}{dz} $$
!
    real, dimension(mz) :: U
    !! Mean intensity array.
    real, dimension(nz) :: dU
    !! First derivative of the mean intensity.
    type (pillar_case) :: p
    !! Pillar object containing the radiative transfer quantitites.
    real :: dz1
    !! Inverse of the grid size along z: \(1/dz\)
!
    ! Calculate dU
    call der(U,dU)
    ! Calculate flux
    p%flux = -1*p%opacity1*dU*dz1

  endsubroutine calc_flux
!====================================================================================
  subroutine get_tridag_coefficients(aa,bb,cc,dd,p,dz,dz2)
    !! Fill in the coefficients of the system of equations in U.
    !! @note
    !! Center and boundary coefficients are calculated differently.
    !! @endnote
!
    real, dimension(nz), intent(inout) :: aa,bb,cc,dd
    !! Tridiag system coefficients
    real, intent(in) :: dz,dz2
    !! Grid size (dz) and its square (dz2) along z.
    type (pillar_case) :: p
    !! Pillar object to update.
!
    call fill_center_coeffs(aa,bb,cc,dd,p,dz2)
    call fill_boundary_coeffs(aa,bb,cc,dd,p,dz,dz2)
!
  endsubroutine get_tridag_coefficients
!====================================================================================
  subroutine fill_center_coeffs(aa,bb,cc,dd,p,dz2)
    !! Populate the center coefficients of the tridiagonal system
    !! for zones i=2 to i=nz-1.
    ! !! i.e.:
    ! !!         _                                          _      _        _
    ! !!       \|    \(a_2\)    \(b_2\)    \(c_2\)       \(0\)       ...     \|  \|   \(d2\)     \|
    ! !!       \|           ...                              \| =\|   ...    \|
    ! !!       \|                  ...                       \| =\|   ...    \|
    ! !!       \|                           ...              \| =\|   ...    \|
    ! !!       \|_   ...     \(0\)    \(a_{nz-1}\)   \(b_{nz-1}\)   \(c_{nz-1}\)  _\|  \|_ \(d_{nz-1}\) _\|
    !! @note
    !! On notation below: `kappa_*` have the same units as `p%opacity`.
    !! @endnote
!
    real, dimension(nz), intent(inout) :: aa,bb,cc,dd
    !! Tridiagonal system coefficients
    real, dimension(nz) :: kappa_m,kappa_p
    !! Opacities at the half points below (`kappa_m`) and above (`kappa_p`).
    real, intent(in) :: dz2
    !! Square of the grid size along z.
    real :: zeta
    !! Intermediate variable for easy reading.
    integer :: iz
    !! Loop variable over z.
    type (pillar_case) :: p
    !! Pillar object.
  
    ! Populate opacities at point i+1/2
    do iz=1, nz-1
      kappa_p(iz) = 0.5 * (p%opacity(iz+1) + p%opacity(iz))
    enddo
     
    ! Populate opacities at point i-1/2
    do iz=2, nz
      kappa_m(iz) = 0.5 * (p%opacity(iz) + p%opacity(iz-1))
    enddo

    ! Populate centers of arrays 
    ! Each coefficient is multiplied by opacity**2 to match the order
    !   of the boundary coefficients.
    do iz=2, nz-1
      aa(iz) = p%opacity2(iz) / kappa_m(iz)
      cc(iz) = p%opacity2(iz) / kappa_p(iz)
      zeta   = dz2 * p%opacity(iz)**3 * (1 - p%albedo(iz))
      bb(iz) = -(aa(iz) + cc(iz) + zeta)
      dd(iz) = -p%source_function(iz) * zeta
    enddo

    if (lfirst) print*, ' Filled center coeffs'
    
  endsubroutine fill_center_coeffs
!====================================================================================
  subroutine fill_boundary_coeffs(aa,bb,cc,dd,p,dz,dz2)
    !! Populate the center coefficients of the tridiagonal system
    !! for zones i=1 and i=nz. 
    !! i.e.:
    !!        \(a_1\) [  \(b_1\)  \(c_1\)    0    ...   ]       =  [   \(d_1\)   ]
    !! and
    !!            [  ...   0   \(a_{nz}\)  \(b_{nz}\)   ] \(c_{nz}\)  =  [  \(d_{nz}\)  ]
!
    real, dimension(nz), intent(inout) :: aa,bb,cc,dd
    !! Tridiagonal system coefficients.
    real, intent(in) :: dz,dz2
    !! Grid size (dz) and its square (dz2) along z.
    real :: zeta
    !! Intermediate variable for easy reading.
!
    type (pillar_case) :: p
    !! Pillar object
  
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
!====================================================================================
endmodule Auxiliary
