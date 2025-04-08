#!/usr/bin/python3

"""
Plot the 2d Flux emerging from the domain.
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
savedir = "./figs_lcf_nodisk/"
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
    os.makedirs(savedir)

# Has columns [snapshot, time (s), Luminosity (erg/s)]
lightcurve = np.loadtxt(Path(datadir, "lightcurve.txt"))
# limit_idx = 41
lightcurve = lightcurve
time_idx = 1
lum_idx = 2

fluxes = []

filelist = glob.glob(os.path.join(datadir, "flux*.npz"))
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
    # if i == 1:
        # break
fluxes = np.asarray(fluxes)


    #   if i == 5:
        # break

###################################################################
def getlevels(array, resolution, log10=True):
    if log10 == True:
        arraymin = np.nanmin(np.log10(array))
        arraymax = np.nanmax(np.log10(array))
        return np.linspace(arraymin, arraymax, resolution)
    else:
        arraymin = np.nanmin(array)
        arraymax = np.nanmax(array)
        return np.linspace(arraymin, arraymax, resolution)
###################################################################

# Set up the colormap to start at the blue instead of black.
cmap = plt.get_cmap("plasma")
shift_fraction = 0.075
new_cmap = mcolors.LinearSegmentedColormap.from_list(
    'shifted_gnuplot2',
    cmap(np.linspace(shift_fraction, 1, 256))
)

fluxLevels = getlevels(fluxes, 256, log10=True)



# [blue, red, purple, black]; faded=0.2 alpha, bold=0.8 alpha
cfaded = [
    (0,0,1,0.2),
    (1,0,0,0.2),
    (0.8901960784313725, 0.4666666666666667, 0.7607843137254902, 0.2),
    (0,0,0,0.2)
]
cbold = [
    (0,0,1,0.8),
    (1,0,0,0.8),
    (0.5390625,0.16796875,0.8828125),
    (0,0,0,0.8)
]

# plot elements: lw=linewidth, ms=markersize, mew=markeredgewidth
lw=2
ms=10
mew=0
marker='o'




for snapshot in range(len(filelist)):
    # snapshot=21

    labels = [# r'$g$ = '+f'{g_t[k]:.2f}',
            # r'$r$ = '+f'{r_t[k]:.2f}',
            # r'$g-r$ = '+f'{gr_t[k]:.2f}',
            r"$L$ = "+f"{lightcurve[snapshot,lum_idx]:.2e}"+r" erg s$^{-1}$"
            ]

    legend_elements = [# Line2D([0], [0], color=cfaded[0], lw=lw, label=labels[0], marker=marker, ms=ms, markerfacecolor=cbold[0], markeredgewidth=mew),
                    # Line2D([0], [0], color=cfaded[1], lw=lw, label=labels[1], marker=marker, ms=ms, markerfacecolor=cbold[1], markeredgewidth=mew),
                    # Line2D([0], [0], color=cfaded[2], lw=lw, label=labels[2], marker=marker, ms=ms, markerfacecolor=cbold[2], markeredgewidth=mew),
                    Line2D([0], [0], color=cfaded[3], lw=lw, label=labels[0], marker=marker, ms=ms, markerfacecolor=cbold[3], markeredgewidth=mew),
                    ]
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12.5,5.3))


    ###################
    # LIGHTCURVE #
    ###################
    # This plots the whole timeseries
    ax[0].plot(lightcurve[:,time_idx]/86400, lightcurve[:,lum_idx],
               color=cfaded[3], lw=lw, marker=marker, ms=ms, fillstyle='none')

    # This plots a single, filled point for the current snapshot
    ax[0].plot(lightcurve[snapshot,time_idx]/86400, lightcurve[snapshot,lum_idx],
               color=cbold[3], marker=marker, ms=ms, markeredgewidth=mew)

    ax[0].set(title=r"Luminosity"+f"\n{waves_A[0]:.0f}"+r"$\leq \lambda (\AA) \leq$"+f"{waves_A[-1]:.0f}",
              xlabel="Days",
              ylabel=r"$log_{10}\ L$ [erg s$^{-1}$]",
              yscale="log")
    # ax[0].set_yscale("log")
    ax[0].legend(handles=legend_elements, title=f"time = {lightcurve[snapshot,time_idx]/86400:.1f} days", title_fontsize=14, fontsize=14)


    ###################
    # FLUX MAP # 
    ###################


    # cmap = "nipy_spectral_r"
    cmap = "gnuplot2"
    flux_cont = ax[1].contourf(x/H_cm, y/H_cm, np.log10(fluxes[snapshot]), 256,
                               cmap=new_cmap, extend='min', levels=fluxLevels,)
                            #    cmap=new_cmap, vmin=11, levels=fluxLevels,)
    ax[1].set(title='Flux',
            aspect='equal',
            xlabel='x/H',
            ylabel='y/H'
            )
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes(position="right", size="5%", pad=0.05)
    plt.colorbar(flux_cont, label=r"$\log_{10}\ F$ [erg cm$^{-2}$ s$^{-1}$]", orientation="vertical", cax=cax, format="%03.2f")

    plt.subplots_adjust(left=0.08, wspace=0.1, right=0.95)
    

    figname = os.path.join(savedir, f"lcf_{snapshot:04}.png")
    print("Saving figure:", figname)
    plt.savefig(figname, format="png", dpi=800)

    plt.close()
    # break
    
    # sys.exit()

    # plt.show()


