!  -*-f90-*-  (for Emacs)    vim:set filetype=fortran:  (for vim)
  type pillar_case
     real, dimension (nz) :: rho,rho1,T,T1,NHII_NHINHII
     real, dimension (nz) :: number_density,inv_number_density,nHI,nHII,ne,ionization_factor
     real, dimension (nz) :: electron_pressure,e_scatter,theta,theta1,lgtheta,lgtheta2
     real, dimension (nz) :: hm_bf_factor,stim_factor
     real, dimension (nz) :: source_function,albedo,opacity,opacity1,flux
  endtype pillar_case
