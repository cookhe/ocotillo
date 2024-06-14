module Grey

  use Common

  implicit none
  private
  
  public :: grey_parameters
  
contains
!******************************************
!  function  
!  endfunction
!******************************************
  subroutine grey_parameters(rho,T,sigma_grey,&
       B_grey,kappa_grey,omega_grey)

    real, intent(in), dimension(nz) :: rho,T
    real, intent(in) :: sigma_grey
    real, dimension(nz) :: alpha_grey
    real, intent(out), dimension(nz) :: B_grey,kappa_grey,omega_grey
    
    B_grey = sigma_sb*T**4
    alpha_grey = 3.68e22 * T**(-3.5) *  rho
    kappa_grey = (alpha_grey+sigma_grey)*rho
    omega_grey = sigma_grey / (alpha_grey + sigma_grey)

  endsubroutine 
!******************************************
endmodule Grey
