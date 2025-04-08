#!/usr/bin/python3

"""
Plot the 1d flux of the supernova.
"""


import sys
import os
from pathlib import Path
import argparse


from scipy.io import FortranFile
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import numpy as np
import glob
import time

plt.rcParams.update({'font.size': 15})


datadir = "./spec_data/"
savedir = "./flux1d/"
unit_length = 1.477e+14 # cm
Nx = 960
Ny = 960
Lx = 96
Ly = 96
x = np.linspace(-Lx/2, Lx/2, Nx) * unit_length
y = np.linspace(-Ly/2, Ly/2, Ny) * unit_length
H_au = 9.871
H_cm = H_au * 1.495979e+13

# Check for and/or create save location.
if not os.path.isdir(savedir):
    print("Create directory:", savedir)
    os.makedirs(savedir)

# Has columns [snapshot, time (s), Luminosity (erg/s)]
lightcurve = np.loadtxt(Path(datadir, "lightcurve.txt"))
# limit_idx = 41
lightcurve = lightcurve
time_idx = 1
lum_idx = 2


filelist = glob.glob(os.path.join(datadir, "flux*.npz"))
fluxes =  []

starttime = time.time()
for i, file in enumerate(filelist):
    if i%5 == 0:
        print("snapshot", i)
    with np.load(file, 'r') as f:
        if i == 0:
            t = f['t']
            waves_A = f['wavesA']
            waves_cm = waves_A / 1e8
      
        fluxes.append(
            np.trapz(f['flux'], x=waves_cm, axis=0)
        )
    # if i == 0:
    # break
fluxes = np.asarray(fluxes)
print("loadtime", time.time()-starttime)

# Loop through each snapshot and plot the flux along
# the row (x,y)=(0,0) at the center to (Lx,0) at the dge
for snapshot in range(fluxes.shape[0]):

    snaptime = lightcurve[snapshot,time_idx]/86400

    label = f"t={snaptime:.1f} Days"
    plt.plot(x[Nx//2:]/H_cm, fluxes[snapshot,Ny//2,Nx//2:],
             label=label, color='indigo', alpha=0.7, lw=2)
    plt.legend()
    plt.xlabel("x/H")
    plt.ylabel(r"$\log_{10}\ F$ [erg cm$^{-2}$ s$^{-1}$]")
    plt.ylim(fluxes.min()/2, fluxes.max()*10)
    plt.yscale('log')
    plt.grid('on')
    plt.tight_layout()
    figname = os.path.join(savedir, f"flux1d_{snapshot:04}.png")
    print("Saving figure:", figname)
    plt.savefig(figname, format="png")
    # plt.show()
    plt.close()

    # break