!  -*-f90-*-  (for Emacs)    vim:set filetype=fortran:  (for vim)
  type pillar_case
  !! Type for holding quantities along a single z column.
     real, dimension (nz) :: rho,rho1,T,T1,NHII_NHINHII
     !! Gas volume density, inverse of ensity, temperature, inverse of temperature,
     !! and ratio of ionized hydrogen (NHII) to total hydrogen (NHI + NHII)
     real, dimension (nz) :: number_density,inv_number_density,nHI,nHII,ne,ionization_factor
     !! Gas number density, inverse number density, number density of neutral hydrogen,
     !! number density of ionized hydrogen, electron number density, base 10 ionization factor
     real, dimension (nz) :: electron_pressure,e_scatter,theta,theta1,lgtheta,lgtheta2
     !! Electron pressure, electron scattering opacity, energy parameter,
     !! inverse energy parameter, log of energy parameter, and square of log energy parameter
     real, dimension (nz) :: hm_bf_factor,stim_factor
     !! Hydrogen minus bound-free ionization factor, and stimulated emission factor
     real, dimension (nz) :: source_function,albedo,opacity,opacity1,opacity2,flux
     !! Source function, albedo, opacity, inverse opacity, opacity sqared, and flux.
  endtype pillar_case
