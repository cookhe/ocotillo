#!/usr/bin/env python3

from datetime import datetime
import numpy as np
from scipy.interpolate import interp1d
import astropy.constants as const

# local imports
import disk
import continuous_opacity as co
import gas_state as gs


def arg():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--density",
        help="Gas density in cgs units (g cm^-3). Default: 1e-9 g cm^-3",
        type=np.float64,
        default=1e-9,
        required=False
    )
    parser.add_argument(
        "-t", "--temperature",
        help="Gas temperature in Kelvin",
        type=np.float64,
        default=1000,
        required=False
    )
    parser.add_argument(
        "-n", "--nlambda",
        help="Integer number of wavelengths to be used. Default: nlambda = 2",
        type=int,
        default=2,
        required=False
    )
    parser.add_argument(
        "--lambdaMin",
        help="Minimum wavelength in Angstroms. Default: lambdaMin = 3,000 A",
        required=False,
        type=np.float64,
        default=3e3
    )
    parser.add_argument(
        "--lambdaMax",
        help="Maximum wavelength in Angstroms. Default: lambdaMax = 10,000 A",
        required=False,
        type=np.float64,
        default=1e4
    )
    parser.add_argument(
        "-v", "--verbose",
        help="If set, print verbose output",
        action="store_true",
    )
    opts = parser.parse_args()

    return opts


def main():
    """
    """

    args = arg()

    c = const.c.cgs.value
    mp = const.m_p.cgs.value
    mp1  = 1/mp
    kb = const.k_B.cgs.value
    sigma_sb = const.sigma_sb.cgs.value


    # define wavelength grid
    lambda_min = args.lambdaMin
    lambda_max = args.lambdaMax
    num_waves = args.nlambda
    waves_angstrom = np.linspace(lambda_min, lambda_max, num_waves)
    waves_cm = waves_angstrom / 1e8 # convert angstroms to centimeters
    nu_Hz = c / waves_cm
    print(f'Calculating for {num_waves} wavelength bins over the range [{lambda_min:.0f}, {lambda_max:.0f}] Angstroms.')

    # Replace isothermal disk with stratified one.
    # tempfloor = 2000
    # T = disk.strat_temp(temp[:,iy,ix], z, floorTemp=tempfloor)
    T = args.temperature # K

    rho = args.density # g cm^-3
    # rhofloor = 1e-15
    # rho = disk.apply_floor(rho, rhofloor)

    # Total number density if all hydrogen
    n = rho * mp1
    # Use Saha Equation to get number densities of HI, HII, and electrons
    nHI, nHII, ne, ionfraction = gs.solve_number_densities(rho, T, n)
                                    
    # Initialize containers
    kappa_H_bf = np.zeros(num_waves)
    kappa_H_ff = np.zeros(num_waves)
    kappa_Hm_bf = np.zeros(num_waves)
    kappa_Hm_ff = np.zeros(num_waves)
    kappa_rad = np.zeros(num_waves)
    absorp_coeff = np.zeros(num_waves)
    source = np.zeros(num_waves)
    omega = np.zeros(num_waves)

    # Criterion for switching between two functions for the free-free absorption
    switch_ionfraction = 1e-2 # i.e. value of (1 - NHII_NHINHII)

    # This ionization_factor is necessary to convert certain opaicites from unts of
    # per neutral hydrogen atom to per hydrogen particle
    ionization_factor = 1/(1 + nHII/nHI)

    # Electron scattering : cm^2 / H
    e_scatter = co.get_electron_scattering_xsec(nHII, n)

    # `theta` has units of 1/eV - used in expressions with 10^(E/kT) when E has units of eV
    theta = 5040./T

    # negative hydrogen (hm) bound-free (bf) factor
    electron_pressure = ne * kb * T
    hm_bf_factor = co.get_hm_bf_factor(electron_pressure, theta)


    # Loop through wavelengths to compute opacities
    for iwave in range(num_waves):
        lfirst = iwave==0
        wave_cm = waves_cm[iwave]
        wave_angstrom = waves_angstrom[iwave]
                        
        ## Compute the source function - Blackbody ##
        source = co.source_function(wave_cm, T)

        # Stimulated emission factor
        stim_factor = co.get_stimulation_factor(waves_angstrom[iwave], theta)

        ## H opacity ##
        # Bound-free neutral hydrogen absorption
        kappa_H_bf[iwave]  = co.xsec_bfHI(np.array([wave_angstrom]), T, m=6)
        kappa_H_bf[iwave] *= stim_factor * ionization_factor

        ## H-minus ion opacity ##
        # Bound-free
        alpha_Hmbf          = co.xsec_bfHminus(np.array([wave_angstrom]))
        kappa_Hm_bf[iwave]  = hm_bf_factor * alpha_Hmbf 
        kappa_Hm_bf[iwave] *= stim_factor * ionization_factor

        # Free-free
        kappa_Hm_ff[iwave]  = co.xsec_ffHminus(ne, np.array([wave_angstrom]), T)
        kappa_Hm_ff[iwave] *= ionization_factor
        
            
        # Total absorption: cm^2 / H
        kappa_rad[iwave] = kappa_H_bf[iwave] + kappa_H_ff[iwave] + kappa_Hm_bf[iwave] + kappa_Hm_ff[iwave]

        # Add electron scattering
        kappa_rad[iwave] += e_scatter

        # Convert to cm^2 / g
        kappa_rad[iwave] /= mp

        # Opacity: 1/cm
        absorp_coeff[iwave] = kappa_rad[iwave] * rho
        if args.verbose and lfirst:
            print()
            print("escatter")
            print("  ",e_scatter)
            print("rho")
            print("  ",rho)#[i1],rho[i2],rho[i3])
            print("absorp_coeff before Hff:")
            print("  ",absorp_coeff[iwave])

        # Free-free absorption
        if 1-ionfraction > switch_ionfraction:
            # Use Gray 2022 function that depends on hydrogen's ionization state.
            kappa_H_ff[iwave]    = co.xsec_ffH(np.array(wave_angstrom), T)
            kappa_H_ff[iwave]   *= stim_factor * ionization_factor # cm^2 / H
            kappa_rad[iwave]    += kappa_H_ff[iwave] / mp # cm^2 / g
            absorp_coeff[iwave] += kappa_rad[iwave] * rho # 1/cm
        else:
            # Use espression for fully ionized gas.
            absorp_coeff[iwave] += co.bremsstrahlung_absorptionCoeff(nu_Hz[iwave], T, ne, nHI+nHII, 1) # 1/cm
            kappa_rad[iwave]    += absorp_coeff[iwave] / rho # cm^2 / g
            kappa_H_ff[iwave]    = kappa_rad[iwave] * mp # cm^2 / H

        # Single scattering albedo: omega = sca / (abs + sca), and match the units again.
        omega[iwave] = (e_scatter) / (kappa_rad[iwave] * mp)

        if args.verbose and lfirst:
            print()
            print("kappa_H_bf")
            print("  ",kappa_H_bf[iwave])
            print("kappa_Hm_bf")
            print("  ",kappa_Hm_bf[iwave])
            print("kappa_Hm_ff")
            print("  ",kappa_Hm_ff[iwave])
            print("kappa_H_ff")
            print( kappa_H_ff[iwave])
            print("K_rad before Hff:")
            print("  ",kappa_rad[iwave])
            print("K_rad after Hff:")
            print("  ",kappa_rad[iwave])
            print("absorp_coeff after Hff")
            print("  ",f"{absorp_coeff[iwave]:14.8e}")

    # The grey case for comparison
    B_grey = sigma_sb*T**4
    alpha_grey = 3.68e22 * T**(-3.5) *  rho
    sigma_grey = 0.4
    kappa_grey = (alpha_grey+sigma_grey)*rho
    absorp_coeff_grey = kappa_grey
    omega_grey = sigma_grey / (alpha_grey + sigma_grey)

    print("")
    print("Gas thermal properties:")
    print(f"Temperature: {T:.2f} K")
    print(f"Density: {rho:.2e} g cm^-3")
    print("")
    print("Continuous absorption coefficients by wavelength:")
    print(f"{'Wavelength (A)':>14s} {'Absorp. Coeff. (1/cm)':>25s}")
    for iwave in range(num_waves):
        print(f"{waves_angstrom[iwave]:14.2f} {absorp_coeff[iwave]:25.15e}")

    print(f"{'Grey':>14s} {absorp_coeff_grey:25.15e}")

if __name__ == "__main__":
    main()      
