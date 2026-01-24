import numpy as np
import astropy.constants as const

#==================================================================

def calc_phiHminus(temps):
    """Phi value for Saha equation (needs to be divided by electron density).
    Parameters: temps : ndarray
                    array of temperatures in Kelvin
    """
    import numpy as np
    
    k        = 1.380649e-16    # [erg K^-1]: Boltzmann
    b        = 4.83e15    # cm^-3 K^-3/2 : phys. constants
    U_Hminus = 1          # partition function for H-
    U_HI     = 2          # partition function for HI
    
    # excitation energy of H-
    chi_Hminus = 0.755 * 1.60218e-12 # convert eV to erg
    
    phi_Hminus = b * temps**(3/2) * (U_HI/U_Hminus) * np.exp(-chi_Hminus/k/temps)
    
    return phi_Hminus

#==================================================================

def fullyionized_ne(n, nelec, abundFrac):
    """Calculate electron number density if all species
    are fully ionized.
    Parameters: n : float
                    total particulate number density
                nelec : array, integer
                    array or list of the total number of 
                    electrons each species contributes to
                    the particle field when fully ionized
                abundFrac : array, float
                    abundance fractions of species
    """
    import numpy as np
    
    ne_nN = 0
    for i, eta in enumerate(abundFrac):
        ne_nN += nelec[i] * eta

    nN = n / (1 + ne_nN)

    ne_max = ne_nN * nN
    return ne_max

#==================================================================

def hydrogen_ionization_fraction(density, temperature):
    """Calculate the fraction of ionized hydrogen HII to HI + HII.

    Args:
        density (float, ndarray): Density in cgs units. Can be a single value or an array.
        temperature (float, ndarray): Gas temperature in cgs units. Can be a single value or an array.
    
    Returns:
        float, ndarray: Fraction of HII abundance to the total abundance of HI + HII.
    """
    mp    = const.u.cgs.value
    me    = const.m_e.cgs.value
    k_eV  = const.k_B.to('eV/K').value
    k_cgs = const.k_B.cgs.value
    h     = const.h.cgs.value

    # calculate the Saha equation (relative fraction of adjacent ions)
    hydrogen_ionization_eV = 13.6 # ionization potential energy for Hydroge in eV
    constants = ((np.sqrt(2*np.pi*me*k_cgs*temperature)/h)**3)#.decompose().cgs
    niine_ni =  constants * np.exp(-hydrogen_ionization_eV/(k_eV*temperature))
    
    C = mp / density * niine_ni

    fraction_HII = (np.sqrt(C**2 + 4*C) - C)/2.

    # Limits due to machine precision mean 4C becomes unresolved next to C**2,
    # so manually set the ionization fraction to 1 for temperatures above 20000.
    # This is appropriate for our domain where the density is always less than
    # ~2e-9 g cm^-3. Would not be appropriate for higher densities.
    # fraction_HII[np.where(temperature>20000)] = 1
    fraction_HII = np.where(temperature>20000, 1, fraction_HII)
    
    return fraction_HII

#==================================================================

def solve_number_densities(density, temperature, number_density):
    """Calculate the number densities of neutral and ionized hydrogen
    and resulting free electrons.
    Parameters
    ----------
    density : float, ndarray
        Gas density.
    temperature : floar, ndarray
        Gas temperature.
    
    Returns
    -------
    tuple of ndarrays
        Number densities of neutral hydrogen (HI), ionized hydrogen (HII),
        free electrons (ne), and ionization fraction (HII / (HI + HII)).
    """

    # Calculate the fraction of HII to total Hydrogen (N = NHI + NHII)
        # i.e., the Saha Equation
    NHII_NHINHII = hydrogen_ionization_fraction(density, temperature)
                
    # Number densities of HII, HI, and electrons
    # APPROXIMATION: all electrons contributed by HII remain free
    nHII = NHII_NHINHII * number_density
    nHI = number_density - nHII
    ne = nHII

    return nHI, nHII, ne, NHII_NHINHII