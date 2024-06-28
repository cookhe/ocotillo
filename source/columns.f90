module Columns

  use Common

  implicit none

  public
  
  integer, parameter :: ncolumns=19
  type column_case
     real, dimension (nz) :: rho,rho1,T,T1,NHII_NHINHII
     real, dimension (nz) :: number_density,inv_number_density,nHI,nHII,ne,ionization_factor
     real, dimension (nz) :: electron_pressure,e_scatter,theta,theta1,lgtheta,lgtheta2
     real, dimension (nz) :: hm_bf_factor,stim_factor
  endtype column_case
  character (len=18), parameter, dimension(ncolumns) :: column_names = (/ &
       'rho               ', &
       'rho1              ', &
       'T                 ', &
       'T1                ', &
       'NHII_NHINHII      ', &
       'number_density    ', &
       'inv_number_density', &
       'nHI               ', &
       'nHII              ', &
       'ne                ', &
       'ionization_factor ', &
       'electron_pressure ', &
       'e_scatter         ', &
       'theta             ', &
       'theta1            ', &
       'lgtheta           ', &
       'lgtheta2          ', &
       'hm_bf_factor      ', &
       'stim_factor       ' &
       /)

endmodule Columns
