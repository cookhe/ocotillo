#!/usr/bin/python3

import sys
import os
from pathlib import Path
import argparse


from scipy.io import FortranFile
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
import numpy as np
import glob
import time

plt.rcParams.update({'font.size': 15})

"""
Calculate the emerging spectral flux and total luminosity from the 
output of an Ocotillo run.
"""

# parse command line arguments
ap = argparse.ArgumentParser()
ap.add_argument("--out",
                 help="Required - Path to directory for saving \
                       concatenated files. Default: current directory.",
                 type=str,
                 default="./",
                 required=False
                )
ap.add_argument("--start",
                 help="Snapshot number (integter) at which to start. Default = 0",
                 type=int,
                 required=False,
                 default=0
                )
ap.add_argument("--stop",
                 help="Snapshot number (integer) at which to stop - exclusive. \
                       Default: start + 1",
                 type=int,
                 required=False
                )
ap.add_argument("--zpov",
                 help="Select the point of view of the observer along the z-axis. \
                  0: the observer views from below. -1: the observer views from above.",
                 type=int,
                 required=True
                )
args = ap.parse_args()
start = args.start
if args.stop:
    stop = args.stop
else:
    stop = start + 1
ns = stop - start

# specify data location and snapshot
datadir0 = "../ocotillo_blue/output/"
datadir1 = "../ocotillo_red/output/"
t = np.loadtxt("../athena/snapshot_times.txt")
# save locations
saveloc = args.out
lcSavePath = Path(saveloc, "lightcurve.txt")
fluxSavePath = Path(saveloc)
# specSavePath = Path(saveloc, "spectral_evolution.txt")

# check if save destination exists, create if not
try:
    os.makedirs(saveloc, exist_ok=False) # won't clobber
    print(f'Directory \'{saveloc}\' created...')
except FileExistsError:
    print(f'Directory \'{saveloc}\' exists...')

# Make luminosity file
if not lcSavePath.is_file():
    with open(lcSavePath, 'a') as f:
        f.write("# snapshot\ttime [s]\tL [erg/s]\n")

# Set the observer's viewing direction
zpov = args.zpov


# set the following parameters to match those used to run
#   the athena and radiative transfer calculations
# Athena:
nprocy = 16    # domain decomposition along y
nprocx = 16   # domain decomposition along x
# nz = 640     # number of grid zones in z
ny = 960      # number of grid zones in y
nx = 960      # number of grid zones in x

# RT:
ng = 3        # number of ghost zones


# timing
# dt = 0.001 # base snapshot output interval from athinput file
# dts = dt*unit_time*snapstep # time elapsed between the selected snapshots (seconds)
# endtime = dts*(ns-1) # final time (assumes starting at 0)
# t = np.arange(0, endtime+dts, dts) # time of each snapshot
unit_time = 44054270.2577661 # s
t = t * unit_time

# set cadence to 0 to turn off or a value indicating the number of 
# wavelengths_angstrom to skip before plotting the next one. 
wl_cadence = 1

# open the z grid data including ghost zones
zfile = FortranFile(os.path.join(datadir0, "zgrid.bin"), "r")
z = zfile.read_record(dtype='d')
zfile.close()
# calculate derived quantities
mz = len(z) # includes ghost zones
nz = mz - 2*ng
dz = z[ng+1]-z[ng]
dx = dz
dy = dz
n0 = ng
n1 = mz - ng
nxloc = nx // nprocx
nyloc = ny // nprocy

######################################################
def read_waves(filepath):
    """Open fortran binary wavelengths."""
    
    wavefile = FortranFile(os.path.join(filepath, "wavegrid.bin"), "r")
    wavelengths_angstrom = wavefile.read_record(dtype='d')
    wavefile.close()
    return wavelengths_angstrom
######################################################
def read_flux_reshape(filepath, nw):
    """Open fortran binary Flux file and reshape data to specifications."""
    
    file = FortranFile(filepath, "r")
    data = file.read_record(dtype='d')
    data = data.reshape(nw,nxloc,nyloc,nz)
    file.close()
    return data
######################################################
def plot_lightcurve(savelocation):

    plt.plot(snapshot_times/86400, luminosity, color="k")
    plt.ylabel(r"log$_{10}\ L$ (erg s$^{-1}$)")
    plt.xlabel(r"Time (days)")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(os.path.join(savelocation, "lightcurve.png"))
    plt.close()
######################################################
def plot_spectrum_evolution(savelocation, cadence=1):

    flux = flux_lam * wavelengths_angstrom
    tdays_map = np.round(snapshot_times/86400/10)*10
    norm = mpl.colors.Normalize(vmin=tdays_map.min(), vmax=tdays_map.max())
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.YlOrBr)
    # cmap.set_array([])
           
    fig, ax = plt.subplots()
    for i in range(0,ns,cadence):
        ax.plot(wavelengths_angstrom, flux[i], c=cmap.to_rgba(tdays_map[i]+1))
    
    cbar = fig.colorbar(cmap, ticks=tdays_map)
    cbar.set_label('Days')
    plt.ylabel(r"$\lambda F_\lambda$ (erg cm$^{-2}$)")
    plt.xlabel(r"Wavelength ($\AA$)")
    plt.legend(title='Days')
    plt.tight_layout()
    plt.savefig(os.path.join(savelocation, "spectral_evolution.png"))
    plt.close()
######################################################

# open the wavelength grid data
wavelengths_angstrom0 = read_waves(datadir0)
wavelengths_angstrom1 = read_waves(datadir1)
wavelengths_angstrom = np.concatenate([wavelengths_angstrom0, wavelengths_angstrom1])
nw0 = len(wavelengths_angstrom0)
nw1 = len(wavelengths_angstrom1)
nw = len(wavelengths_angstrom)
# get frist and last wavelength
w0s = wavelengths_angstrom0[0]
w0e = wavelengths_angstrom0[-1]
w1s = wavelengths_angstrom1[0]
w1e = wavelengths_angstrom1[-1]
ws = wavelengths_angstrom[0]
we = wavelengths_angstrom[-1]
wavelengths_cm = 1e-8 * wavelengths_angstrom


# set up shapes of containers
# flux_shape = (ns, nw, ny, nx)
flux_shape = (nw, ny, nx)
pillar_shape    = (nw, nxloc, nyloc, nz)
pillar_shape_gz = (nw, nxloc, nyloc, mz)
count    = int(np.prod(pillar_shape))
count_gz = int(np.prod(pillar_shape_gz))

# global grids
# wavegrid = np.linspace(w0,w1,nw)

# lines for plot legend
custom_lines = [Line2D([0], [0], color='k', ls='-'),
                Line2D([0], [0], color='k', ls='--')]


# TO-DO: improve by only opening the column slice
# associated with the chosen xloc,yloc coordinates.

start_time = time.time()

# do the loop per frequency as well once the file is open.
# perform the calculations then add the result to the U.
# U array could maybe just be a 2d array in XY. Values 
# calculated from the processor columns are simply added
# to the appropriate location after 
for snapnum  in range(start, stop):
    snapshot_start_time = time.time() # timing for each snapshot loop
    snapshot = f"{snapnum:04}" # number padded with zeros to fill 4 digits
    print(f"Snapshot: {snapshot}")
    flux = np.zeros(flux_shape) # V
    print(" Starting with x-y loop:")
    for ixp in range(nprocx):
        print(f" x proc: {ixp}")
        for iyp in range(nprocy):
            print(f"  y proc: {iyp}")

            file0 = os.path.join(datadir0, f"procx{ixp}_procy{iyp}/flux_{snapshot}.bin")
            file1 = os.path.join(datadir1, f"procx{ixp}_procy{iyp}/flux_{snapshot}.bin")
    
            # Select the data at the boundary closest to the observer's POV.
            V0 = read_flux_reshape(file0, nw0)[:,:,:,zpov]
            V1 = read_flux_reshape(file1, nw1)[:,:,:,zpov]

            # need the indicies for the columns in context of the global arrays
            xpstart = ixp*nxloc
            xpend   = xpstart + nxloc
            ypstart = iyp*nyloc
            ypend   = ypstart + nyloc
            pillar_loc0 = np.s_[  0:nw0, xpstart:xpend, ypstart:ypend]
            pillar_loc1 = np.s_[nw0:   , xpstart:xpend, ypstart:ypend]

            # the emerging flux is the last point in Z
            flux[pillar_loc0] = V0 * 2 * np.pi * np.repeat(wavelengths_cm[:nw0],nxloc*nyloc).reshape(nw0,nyloc,nxloc)
            flux[pillar_loc1] = V1 * 2 * np.pi * np.repeat(wavelengths_cm[nw0:],nxloc*nyloc).reshape(nw1,nyloc,nxloc)

        
            # if iyp == 0:
            #     break
            """end of nyproc loop"""
        # if ixp == 0:
        #     break
        """end of nxproc loop"""

    # save quantities to file with keywords
    fluxSaveName = Path(fluxSavePath, f"flux_{snapshot}.npz")
    print(f'Saving compressed Flux() file: {fluxSaveName}')
    np.savez_compressed(fluxSaveName,
                        flux=flux,
                        t=t[snapnum],
                        wavesA=wavelengths_angstrom
                        )

    # Do flux and luminosity calculations on the full snapshot
    flux_sum = np.sum(flux, axis=(1,2)) # sum over xy
    flux_lam = flux_sum / nx / ny # return to units of erg/s/cm^2/A
    luminosity = np.trapz(flux_sum * dz * dz, x=wavelengths_cm, axis=-1)
    
    # Add a line with this snapshot's luminosity to the lightcurve file
    with open(lcSavePath, 'a') as f:
        print(f"Writing snapshot luminosity to {lcSavePath.as_posix()}")
        f.write(f"{snapshot}\t{t[snapnum]:e}\t{luminosity:e}\n")

    print('Snapshot time', time.time() - snapshot_start_time)
    # break
    """end of snapshot loop"""

print('Total execution time', time.time() - start_time)

# There may be more 
# snapshot_times = t[start:stop]

# lc_header = f"{'time (s)':>8s}{'erg/s':>8s}"
# np.savetxt(lcSavePath, np.c_[snapshot_times,luminosity], header=lc_header, fmt="%.8e")

# spec_header = f"{'time (s)':>12s}   {'erg/s/cm^2/A'}"
# np.savetxt(specSavePath, np.c_[snapshot_times,flux_lam], header=spec_header, fmt="%.8e")

# print some diagnostics
# print("")
# print("flux\n",flux[0,:,0,0])
# print("flux_sum\n",flux_sum[0,:])
# print("nx, ny\n", nx, ny)
# print("flux_lam\n", flux_lam[0,:])
# print("median flux_lam\n", np.median(flux_lam,axis=0))
# print("dz/1e13, (dz*dz)/1e26\n",dz/1e13,(dz*dz)/1e26)
# print("wavelengths_cm\n",wavelengths_cm)
# print("luminosity\n", luminosity[0])


# plot_lightcurve(saveloc)
# plot_spectrum_evolution(saveloc, cadence=wl_cadence)


##########################################################################
# Code to open mean intensity (U) and opacity(k) data

# Ufile_path = f"{datadir}procx{ixp}_procy{iyp}/"+ \
#     f"mean_intensity_{snapshot}.bin"
# print("opening file: ", Ufile_path)
# Ufile = FortranFile(Ufile_path, "r")
# U = Ufile.read_record(dtype='d')
# U = U.reshape(nw,nxloc,nyloc,mz)
# Ufile.close()

# kfile_path = f"{datadir}procx{ixp}_procy{iyp}/"+ \
#     f"absorption_coefficients_{snapshot}.bin"
# print("opening file: ", kfile_path)
# kfile= FortranFile(kfile_path, "r")
# k = kfile.read_record(dtype='d')
# k = k.reshape(nw,nxloc,nyloc,nz)
# kfile.close()
