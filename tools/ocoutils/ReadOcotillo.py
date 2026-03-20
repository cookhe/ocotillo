#!/usr/bin/python3

import os
import numpy as np
from scipy.io import FortranFile
import h5py

def openbinary(filepath, verbose=False):
    try:
        file = open(filepath, "rb")
        if verbose:
            print(f"opened {filepath}")
    except:
        print("data must be <binary_dump>")
        print(f"check {filepath}")
        raise SystemExit

    # get the end of file trigger. Then, return to beginning
    file.seek(0,2)
    eof = file.tell()
    file.seek(0,0)

    return file, eof
def der6(f):
    """Compute the sixth-order derivative of the array f.

    Parameters
    ----------
    f : array
        Some array of numbers over which to calculate the derivative.

    Returns
    -------
    array
        The sixth-order derivative of f.
    """
    df=np.zeros(len(f))
    for i in range(3,len(f)-4):
        df[i] = 1/60*(f[i+3]-f[i-3]) - 3./20*(f[i+2]-f[i-2]) + 3./4*(f[i+1]-f[i-1])
    return df


def read_waves(filepath):
    """Open fortran binary wavelengths.

    Parameters
    ----------
    filepath : str
        Path to directory containing wavegrid.bin

    Returns
    -------
    wavelengths_angstrom : ndarray
        Wavelength grid in Angstroms
    """
    wavefile = FortranFile(os.path.join(filepath, "wavegrid.bin"), "r")
    wavelengths_angstrom = wavefile.read_record(dtype='d')
    wavefile.close()
    return wavelengths_angstrom

def read_zgrid(filepath):
    """Open fortran binary zgrid. Here, z is along calculated direction,
    so it may not reflect the z-dimension in the hydrodynamic grid.

    Parameters
    ----------
    filepath : str
        Path to directory containing wavegrid.bin

    Returns
    -------
    z : ndarray
        Grid points along direction used for ocotillo calculation
    """
    # Read zgrid to get dz
    zfile = FortranFile(os.path.join(rtdatadir, "zgrid.bin"), "r")
    z = zfile.read_record(dtype='d')
    zfile.close()
    return z

def read_flux_binary_reshape(filepath, nw, nxloc, nyloc, nz):
    """Open fortran binary Flux file and reshape data to specifications.

    Parameters
    ----------
    filepath : str
        Full path to flux binary file
    nw : int
        Number of wavelengths
    nxloc : int
        Local grid size in x
    nyloc : int
        Local grid size in y
    nz : int
        Grid size in z (no ghost zones)

    Returns
    -------
    data : ndarray, shape (nw, nxloc, nyloc, nz)
        Flux data
    """
    file = FortranFile(filepath, "r")
    data = file.read_record(dtype='d')
    data = data.reshape(nw, nxloc, nyloc, nz)
    file.close()
    return data

def read_mean_intensity_binary_reshape(filepath, nw, nxloc, nyloc, mz):
    """Open fortran binary mean intensity file and reshape.

    Parameters
    ----------
    filepath : str
        Full path to mean_intensity binary file
    nw : int
        Number of wavelengths
    nxloc : int
        Local grid size in x
    nyloc : int
        Local grid size in y
    mz : int
        Grid size in z including ghost zones

    Returns
    -------
    data : ndarray, shape (nw, nxloc, nyloc, mz)
        Mean intensity data (includes ghost zones)
    """
    file = FortranFile(filepath, "r")
    data = file.read_record(dtype='d')
    data = data.reshape(nw, nxloc, nyloc, mz)
    file.close()
    return data

def read_absorption_coefficients_binary_reshape(filepath, nw, nxloc, nyloc, nz):
    """Open fortran binary absorption coefficient file and reshape.

    Parameters
    ----------
    filepath : str
        Full path to absorption_coefficients binary file
    nw : int
        Number of wavelengths
    nxloc : int
        Local grid size in x
    nyloc : int
        Local grid size in y
    nz : int
        Grid size in z (no ghost zones)

    Returns
    -------
    data : ndarray, shape (nw, nxloc, nyloc, nz)
        Absorption coefficient data
    """
    file = FortranFile(filepath, "r")
    data = file.read_record(dtype='d')
    data = data.reshape(nw, nxloc, nyloc, nz)
    file.close()
    return data

def get_column_data(datadir, snapshot, x_choice, y_choice,
                    nprocx, nprocy, nxloc, nyloc, nz, nw, ng=3,
                    data_type='absorption_coefficients'):
    """Get mean intensity or absorption coefficient as function of z and wavelength
    at a specific xy coordinate.

    Parameters
    ----------
    datadir : str
        Base directory containing processor subdirectories
    snapshot : str or int
        Snapshot number (will be formatted as 4-digit string)
    x_global : int
        Global x index (0 to nx-1)
    y_global : int
        Global y index (0 to ny-1)
    nprocx : int
        Number of processors in x direction
    nprocy : int
        Number of processors in y direction
    nxloc : int
        Local grid size in x per processor
    nyloc : int
        Local grid size in y per processor
    nz : int
        Grid size in z (no ghost zones)
    nw : int
        Number of wavelengths
    ng : int, optional
        Number of ghost zones (default: 3)
    data_type : str, optional
        Either 'absorption_coefficients' or 'mean_intensity' (default: 'absorption_coefficients')

    Returns
    -------
    data : ndarray, shape (nw, nz) or (nw, mz)
        Column of data at the specified xy coordinate
        For mean_intensity, includes ghost zones (mz = nz + 2*ng)
        For absorption_coefficients, no ghost zones (nz)
    """
    # Format snapshot as 4-digit string
    if isinstance(snapshot, int):
        snapshot = f'{snapshot:04}'

    # Determine which processor owns this coordinate
    ixp = x_choice // nxloc
    iyp = y_choice // nyloc

    if ixp > nprocx-1 or iyp > nprocy-1:
        e = f"Coordinates (x_choice, y_choice)=({x_choice}, {y_choice}) outside " \
            f"domain with dimensions (nx, ny)=({nxloc*nprocx}, {nyloc*nprocy})"
        raise ValueError(e)
    if ixp < 0 or iyp < 0:
        e = f"Coordinates must be positive but are " \
            f"(x_choice, y_choice)=({x_choice}, {y_choice})"
        raise ValueError(e)

    # Local indices within that processor
    x_local = x_choice - ixp * nxloc
    y_local = y_choice - iyp * nyloc


    # Construct file path
    if data_type == 'mean_intensity':
        filename = f'mean_intensity_{snapshot}.bin'
        mz = nz + 2 * ng
        filepath = os.path.join(datadir, f'procx{ixp}_procy{iyp}', filename)
        data = read_mean_intensity_binary_reshape(filepath, nw, nxloc, nyloc, mz)
        # Extract the column
        column = data[:, x_local, y_local, :]  # has shape (nw, mz)
    elif data_type == 'absorption_coefficients':
        filename = f'absorption_coefficients_{snapshot}.bin'
        filepath = os.path.join(datadir, f'procx{ixp}_procy{iyp}', filename)
        data = read_absorption_coefficients_binary_reshape(filepath, nw, nxloc, nyloc, nz)
        # Extract the column
        column = data[:, x_local, y_local, :]  # has shape (nw, nz)
    elif data_type == 'flux':
        filename = f'flux_{snapshot}.bin'
        filepath = os.path.join(datadir, f'procx{ixp}_procy{iyp}', filename)
        data = read_flux_binary_reshape(filepath, nw, nxloc, nyloc, nz)
        column = data[:, x_local, y_local, :]
    else:
        raise ValueError(f"Unknown data_type: {data_type}. Use 'mean_intensity', 'absorption_coefficients', or 'flux'")

    return column

def calculate_optical_depth(absorption_coeff, density, dz):
    """Calculate optical depth from absorption coefficients.

    Parameters
    ----------
    absorption_coeff : ndarray, shape (nw, nz)
        Absorption coefficients along a column
    dz : float
        Grid spacing in z direction

    Returns
    -------
    tau : ndarray, shape (nw, nz)
        Cumulative optical depth from z=0
    """
    tau = np.cumsum(absorption_coeff * density * dz, axis=1)
    return tau

def find_photosphere(tau, target_tau=1.0):
    """Find z-index where optical depth reaches a target value.

    Parameters
    ----------
    tau : ndarray, shape (nw, nz)
        Optical depth array
    target_tau : float, optional
        Target optical depth (default: 1.0 for photosphere)

    Returns
    -------
    z_indices : ndarray, shape (nw,)
        Z-index where tau crosses target_tau for each wavelength
        Returns -1 if target_tau never reached
    """
    nw, nz = tau.shape
    z_indices = np.zeros(nw, dtype=int)

    for iw in range(nw):
        # Find first index where tau >= target_tau
        idx = np.where(tau[iw, :] >= target_tau)[0]
        if len(idx) > 0:
            z_indices[iw] = idx[0]
        else:
            z_indices[iw] = -1  # Never reached target

    return z_indices

def save_hdf5(filename):


                