import sys
import numpy as np
import astropy.constants as const

c = const.c.cgs.value
h = const.h.cgs.value
kb = const.k_B.cgs.value
e = const.e.esu.value
me = const.m_e.cgs.value


#==================================================================

def gaunt(n, waves):
    """Calculate Gaunt Factor
    Parameters: n : int
                    principle quantum number
                waves : ndarray
                    wavelengths in Angstroms
    """
    
    c = 2.9979246e+10 # cm s^-1 - lightspeed
    h = 6.6260755e-27 # erg s   - Planck constant
    R = 2.1798741e-11 # erg     - Rydberg energy
    
    # convert Angstroms to centimeters
    waves = waves*1e-8
    
    # reference wavelength lambda_n
    lam_n = n**2 * h*c/R
    
    # gaunt will be non-zero only where wavelengths are shorter than
    #   lam_n (i.e. the ionizing wavelength)
    lam_mask = waves <= lam_n
    
    # get epsilon_n value
    eps_n = (waves[lam_mask]/lam_n) - 1
    
    # Gaunt factor. Below experession is for lamda <= lambda_n; zero otherwise
    gaunt_factor = np.zeros(len(waves))
    gaunt_factor[lam_mask] = 1 - (121./700.) * (1-eps_n**2) / (n*(1+eps_n**2))**(2/3)
    
    return gaunt_factor

#==================================================================

def xsec_bfHI(waves, temp, m=6):
    """Cross section of bound-free hydrogen. Sums over the first 
    1 to m-1 excitation states, where m is the principal quantum
    number and uses the Unsold integral approximation for levels
    m to infinity. Returns units of cm^2 per neutral hydrogen atom.

    Parameters
    ----------
    waves : float
        Wavelength in Angstroms
    temp : float
        Temperature in Kelvin
    m : int, optional
        Principal quantum number at which to begin the Unsold
        integral approximation, by default 6

    Returns
    -------
    float
        Cross section of bound-free hydrogen in units of cm^2 per
        neutral hydrogen atom
    """
    
    A = 1.0449e-26    # cm^2 A^-3 - fundamental constants
    k = 1.380658e-16  # erg K^-1  - Boltzmann constant
    R = 2.1798741e-11 # erg       - Rydberg energy
           
    # sum for the first m-1 excitation states
    sm = 0.
    for n in range(1,m):
        chi_n = R * (1 - 1/n**2)
        g     = gaunt(n, waves)
        sm   += g/n**3 * np.exp(-chi_n/k/temp)
    
    # excitation energy of state level m
    chi_m = R * (1 - 1/m**2)
    
    # Unsold approximation integral
    C = k*temp/2/R * ( np.exp(-chi_m/k/temp) - np.exp(-R/k/temp) )
    
    alpha = A * waves**3 * (C + sm)
    
    return alpha

#==================================================================

def xsec_ffH(waves, temp):
    """Hydrogen free-free absorption coefficient.
    
    Units cm^2 per neutral hydrogen atom.
    """
    
    k = 1.380658e-16   # erg K^-1  - Boltzmann constant
    R = 2.1798741e-11  # erg       - Rydberg energy
    c = 2.9979246e+10  # cm s^-1   - lightspeed
    h = 6.6260755e-27  # erg s     - Planck constant
    e  = 4.8032068e-10 # [esu]: charge
    me = 9.1093897e-28 # [g]: electron mass
           
    # absorption per electron
    alpha0 = 1.0443e-26
    Rcm = 2 * np.pi**2 * me * e**4 / h**3 / c
    Rangstrom = 1.0968e-3
    # gaunt factor
    theta = 5040/temp
    chi = 1.2398e4 / waves
    log10e = np.log10(np.e)

    g_ff = 1 + 0.3456 / (waves * Rangstrom)**(1/3) * (log10e/theta/chi + 0.5)
    
    I = h*c*Rcm # erg
    I *= 6.2419e+11 # convert to eV
    energyfactor = log10e/2/theta/I * 10**(-theta*I)
    
    return alpha0 * waves**3 * g_ff * energyfactor
    
#==================================================================

def xsec_bfHminus(waves):
    """H- ion bound-free continuum opacity.
    Calculates the sum of the fit function
    
        alpha = 1e-18 * sum_m=0^6(a_m * wave^m)
    
    Parameters: waves : 1D array
                    array of wavelengths. Single values
                    must be contained in an array object.
    
    """
        
    # fit coefficients
    a = np.array([
        +1.99654,
        -1.18267e-6,
        +2.62423e-7,
        -4.40524e-11,
        +3.23992e-15,
        -1.39568e-19,
        +2.78701e-24
    ])
    
    s = 0 * waves
    m = np.arange(0, len(a), 1)
    
    # loop through each wavelength and perform the summation
    for ind, w in enumerate(waves):
        s[ind] = np.sum( a * w**m )

    alpha = 1e-17 * s
    
    # remove unphysical values
    alpha[waves>16000] = 0
        
    return alpha.T

#==================================================================

def xsec_e(ne):
    """Thomson electron scattering cross section per electron"""
    
    c  = 2.99792458e+10 # [cm]: speed of light
    e  = 4.8032068e-10  # [esu]: charge
    me = 9.1093897e-28  # [g]: electron mass
    
    re = e**2 / me / c**2
    
    return 8.*np.pi/3. * ne * re**2

#==================================================================

def xsec_ffHminus(ne, waves, temp):
    """H- ion free-free interaction cross section per HI per unit 
    electron pressure.
    """
        
    k = 1.380658e-16  # erg K^-1  - Boltzmann constant
    # fit coefficients
    b = np.array([
        [-2.276300,
         -1.685000,
         +0.766610,
         -0.053356,
         +0.000000],
        [+15.28270,
         -9.284600,
         +1.993810,
         -0.142631,
         +0.000000],
        [-197.789,
         +190.266,
         -67.9775,
         +10.6913,
         -0.62515]
    ])
    
    # integers used for exponents
    ks = np.arange(0, len(b[0]), 1)
    
    # make wavelength array for each k
    waves = np.stack((waves, waves, waves, waves, waves), axis=-1)
    # print(waves.shape)
    # 
    f = np.zeros((len(b), len(waves)))
    for i in range(len(b)):
        f[i] = np.sum(np.log10(waves)**ks * b[i], axis=1)
        
    theta = 5040./temp
    
    flam = f[0] + f[1]*np.log10(theta) + f[2]*np.log10(theta)**2 - 26
    
    alpha = ne * k * temp * 10**(flam)
    
    return alpha

def xsec_hyd(waves, T, ne, m):
    '''Calculate the following cross sections:
    HI bound-free
    H- bound-free
    H- free-free
    '''
    
    # HI bound-free alpha
    a_bfHI = xsec_bfHI(waves, T, m)
    
    # H- bound-free alpha
    a_bfHminus = xsec_bfHminus(waves)
    phi_Hminus = calc_phiHminus(T)
    a_bfHminus = (a_bfHminus * ne/phi_Hminus).T
    a_bfHminus[a_bfHminus<0] = 0 # set negative values to zero
    
    # H- free-free alpha
    a_ffHminus = xsec_ffHminus(ne, waves, T)
    
    return np.array([a_bfHI, a_bfHminus, a_ffHminus])

#==================================================================

def selfcorrection(waves, absorber, T):
    '''Self-emission correction for opacity.'''
    
    
    c = 2.9979246e+10 # cm s^-1   - lightspeed
    h = 6.6260755e-27 # erg s     - Planck constant
    k = 1.380658e-16  # erg K^-1  - Boltzmann constant

    # self-emission correction
    K = 1 - np.exp(-h*c/waves/k/T)
        
    return K*absorber

#==================================================================

def gaunt_factor_hertz(nu, temperature):
    """Quantum mechanical correction of order unity."""
    
    return np.log(np.exp(5.960 - np.sqrt(3.)/np.pi * np.log(nu/1e9 * (temperature/1e4)**(-1.5))) + np.exp(1.))

#==================================================================

def bremsstrahlung_absorptionCoeff(frequency, temperature, nelectrons, nprotons, zprotons):

    ln_const = np.log(4/3) + 6*np.log(e) - np.log(me*h*c) + 0.5*np.log(2*np.pi) - 0.5*np.log(3*me*kb)

    gaunt = gaunt_factor_hertz(frequency, temperature)

    ln_physQuant = np.log(nelectrons*nprotons) + 2*np.log(zprotons) - 0.5*np.log(temperature) - 3*np.log(frequency) + np.log((1 - np.exp(-h*frequency/kb/temperature)))

    kappaRho = np.exp(ln_const) * np.exp(ln_physQuant) * gaunt

    return kappaRho

#==================================================================

def get_hm_bf_factor(electron_pressure, theta):
    """Calculate H- bound-free factor used in opacity calculations.
    Parameters
    ----------
    electron_pressure : float, ndarray
        Electron pressure in cgs units.
    theta : float, ndarray
        5040 / temperature in Kelvin.
    
    Returns
    -------
    float, ndarray
        H- bound-free factor.
    """
    return 4.158e-10 * electron_pressure * theta**(5./2) * 10**(0.754 * theta)

#==================================================================

def source_function(wave_cm, T):
    """Calculate the blackbody source function at a given wavelength and temperature.

    Parameters
    ----------
    wave_cm : float, ndarray
        Wavelength in centimeters.
    T : float, ndarray
        Temperature in Kelvin.
    """

    # Must protect the exponential in the source function calculation from overflow.
    # 2**1024 is the maximum size for floats.
    overflow_limit = int(np.floor(np.log10(sys.float_info.max)))

    damping_factor = h*c/wave_cm/kb/T
                
    # Set exceeding values to just below the limit.
    damping_factor = np.where(damping_factor > overflow_limit, overflow_limit, damping_factor)
    return 2*h*c**2/wave_cm**5 * 1/(np.exp(damping_factor)-1)

#==================================================================

def get_electron_scattering_xsec(nHII, n):
    """Electron scattering in cm^2 / H. 

    Parameters
    ----------
    nHII : float, ndarray
        Number density of ionized hydrogen (HII).
    n : float, ndarray
        Total number density.

    Returns
    -------
    float, ndarray
        Electron scattering cross section in cm^2 / H.
    """

    alpha_e = 0.6648e-24 # coefficient
    return alpha_e * (nHII / n)

#==================================================================

def get_stimulation_factor(wave_angstrom, theta):
    """The stimulated emission factor.

    Parameters
    ----------
    wave_angstrom : float, ndarray
        Wavelength in Angstroms.
    theta : float, ndarray
        5040 / temperature in Kelvin.

    Returns
    -------
    float, ndarray
        The stimulated emission factor.
    """
    chi = 1.2398e4 / wave_angstrom
    return 1-10**(-chi*theta)