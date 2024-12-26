#!/usr/bin/python3

import sys

import numpy as np

# set the following parameters to match those used to run
#   the athena and radiative transfer calculations
# Athena:
nprocx = 2#8
nprocy = 1#8
nz = 36#640
ny = 34#480
nx = 32#480
# radiative transfer:
w0 = 3000 # angstroms
w1 = 5000
nw = 5        # number of wavelengths
ng = 3        # number of ghost zones

# specify data
datadir = "/home1/07139/hecook/ocotillo/samples/multiInTest/output/"
snapshot = "0000"


"""
For each snapshot:
    1. create a global 3d grid to hold final calculations: nw, ny, nx
    2. for each processor column, create a local 3d grid: 
        a. mean intensity (U): mz, ny/nprocy, nx/nprocx (mz includes ghost zones)
        b. absorption_coefficient (kappa): nz, ny/nprocy, nx/nprocx
    3. for each wavelength, calculate:
        a. flux (V)
        b. optionally: upward intensity (Ip)
        c. optionally: downward intensity (Im)
    4. save to the respective x,y location in the global grid

"""



mz = nz + 2*ng
nxloc = nx//nprocx
nyloc = ny//nprocy
global_shape = (nw, ny, nx)
pillar_shape    = (nz, nyloc, nxloc, nw)
pillar_shape_gz = (mz, nyloc, nxloc, nw)
count    = int(np.prod(pillar_shape))
count_gz = int(np.prod(pillar_shape_gz))

# get z grid information
zgrid = np.loadtxt(f"{datadir}/zgrid.txt")[:,1]
dz = zgrid[1]-zgrid[0]

# global grids
wavegrid = np.linspace(w0,w1,nw)
flux = np.zeros(global_shape) # V

# pillar grids
# U = np.zeros(pillar_shape_gz)
# abs_coef = np.zeros(pillar_shape)



######################################################
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
######################################################
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
######################################################


# do the loop per frequency as well once the file is open.
# perform the calculations then add the result to the U.
# U array could maybe just be a 2d array in XY. Values 
# calculated from the processor columns are simply added
# to the appropriate location after 
for ixp in range(nprocx):
    for iyp in range(nprocy):
        abs_file_path = f"{datadir}procx{ixp}_procy{iyp}/"+ \
            f"absorption_coefficients_{snapshot}.bin"
        U_file_path = f"{datadir}procx{ixp}_procy{iyp}/"+ \
            f"mean_intensity_{snapshot}.bin"

        abs_file, eof_abs = openbinary(abs_file_path)
        U_file, eof_U = openbinary(U_file_path)
        
        abs_coef = np.fromfile(abs_file, dtype=np.float64, count=count).reshape(pillar_shape)
        U = np.fromfile(U_file, dtype=np.float64, count=count_gz).reshape(pillar_shape_gz)
        # print(abs_coef[:,0,0,0])
        print(U[:,0,0,0])
        # for j in range(nyloc):
        #     for i in range(nxloc):
        #         print(U[:,j,i,0])
                
        # print("min/max(abs_coef)", np.min(abs_coef), np.max(abs_coef))
        # print("min/max(U)", np.nanmin(U), np.nanmax(U))
        # # wavelength loop
        # for iw, iwave in enumerate(wavegrid):

        # # Make variables for global index locations to clean this up
            
            
            # if ixp == 0 and iyp == 0:


            # for 
            # V[:,iwave] = -1/abs_coef * der6(U) / dz
        #     flux[iw, iy*nyloc:(iy+1)*nyloc, ix*nxloc:(ix+1)*nxloc]
            
        #     # reset arrays to zero for next wavelength
        #     U.fill(0)
        #     abs_coef.fill(0)
            # break

        abs_file.close()
        U_file.close()
        break
    break
        # u_file = f"{datadir}procx{ix}_procx{iy}/mean_intensity_{snapshot}.bin"
                