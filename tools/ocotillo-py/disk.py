import numpy as np

def strat_temp(temperature, z, eps=0.2, isoTemp=64857.04298198191, sigma=2.54e14, midTemp=61341, floorTemp=3):
    """Calculate a gausian disk profile and add it to the regions of the provided
    temperature which are not the isothermal temperature

    eps  = threshold percentage above isoTemp. Regions below this value will be rewritten with
           gaussian stratified profile
    isoTemp = 64857.0 K, isothermal temperature of initial conditions
    sigma = 2.54e14 # matches iterativeTemperature-HC.ipynb matching gaussian
    midTemp = 61341 midplane temperature
    floorTemp = temperature floor, default = 3 K

    T (< (1+eps) * isoTemp) = gaussian profile
    """
        
    # gaussian profile based on midplane temperature
    T_gauss = midTemp*np.exp(-0.5*z**2/sigma**2)

    dims = temperature.shape # dimensions of array (assumes 3D)
    # sys.exit()
    # fill in z-columns of 3d array with 1d gaussian profile
    temperature_strat = np.zeros((dims))
    if len(dims) == 1:
        temperature_strat = T_gauss
    elif len(dims) == 2:
        for ix in range(dims[1]):
            temperature_strat[:,ix] = T_gauss
    else:
        for jy in range(dims[1]):
            for jx in range(dims[2]):
                temperature_strat[:,jy,jx] = T_gauss

    # use stratified temp when Athena temp less than 1+eps*isoTemp, use Athena temp when not.
    compositeTemp = (temperature_strat * (temperature <  (1+eps)*isoTemp) +
                           temperature * (temperature >= (1+eps)*isoTemp))
    
    compositeTemp = apply_floor(compositeTemp, floorTemp)

    return compositeTemp

#==================================================================

def apply_floor(array, floor):
    """Apply a floor value to an array

    Parameters
    ----------
    array : ndarray
        Array or scalar containing data to apply the floor.
    floor : ndarray
        Value to set the floor at.

    Returns
    -------
    ndarray
        Array with floor applied.
    """

    return np.where(array < floor, floor, array)
