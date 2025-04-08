#!/usr/bin/python3

"""
Plot the color of the radiation emerging from the domain.
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
savedir = "./figs_color/"
unit_length = 1.477e+14 # cm
Nx = 960
Ny = 960
Lx = 96
Ly = 96
x = np.linspace(-Lx/2, Lx/2, Nx) * unit_length
y = np.linspace(-Ly/2, Ly/2, Ny) * unit_length
H_au = 9.871
H_cm = H_au * 1.495979e+13
Nw=50

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
fluxes = []

for i, file in enumerate(filelist):
    if i%5 == 0:
        print("snapshot", i)
    with np.load(file, 'r') as f:
        if i == 0:
            t = f['t']
            waves_A = f['wavesA']
            waves_cm = waves_A / 1e8
      
        fluxes.append(f['flux'])
    # if i == 1:
    # break
fluxes = np.asarray(fluxes)

# Get the SED and total spatial fluxes
flux_lam = np.sum(fluxes, axis=(2,3)) / Nx / Ny # sum over xy

# Central wavelength of the filter (by eye)
g_center = 4800
r_center = 6300

g_index = np.argmin(np.abs(waves_A-g_center))
r_index = np.argmin(np.abs(waves_A-r_center))


g_flux = flux_lam[:,g_index]
r_flux = flux_lam[:,r_index]

color = -2.5*np.log10(g_flux/r_flux)
deltacolor = color - color[0]

# sys.exit()
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

# plt.(nrows=1, ncols=3, figsize=(18.75,5.3), layout="constrained")
# This plots the whole timeseries
plt.plot(lightcurve[:,time_idx]/86400, color,
        color='k', lw=lw, marker=marker, ms=ms, fillstyle='none')

plt.ylabel(r'$-2.5 \log_{10} \left( \frac{F(4800 \AA)}{F(6300 \AA)} \right)$')
plt.xlabel('Days')
plt.tight_layout()

figname = os.path.join(savedir, f"color.png")
print("Saving figure:", figname)
plt.savefig(figname, format="png")
plt.close()


plt.plot(lightcurve[:,time_idx]/86400, deltacolor,
        color='k', lw=lw, marker=marker, ms=ms, fillstyle='none')

plt.ylabel(r'$(g-r)-(g-r)_0$')
plt.xlabel('Days')
plt.tight_layout()

figname = os.path.join(savedir, f"deltacolor.png")
print("Saving figure:", figname)
plt.savefig(figname, format="png")
plt.close()



# for snapshot in range(len(filelist)):

#     labels = [# r'$g$ = '+f'{g_t[k]:.2f}',
#             # r'$r$ = '+f'{r_t[k]:.2f}',
#             # r'$g-r$ = '+f'{gr_t[k]:.2f}',
#             r"$L$ = "+f"{color[snapshot]:.2e}"+r" erg s$^{-1}$"
#             ]

#     legend_elements = [# Line2D([0], [0], color=cfaded[0], lw=lw, label=labels[0], marker=marker, ms=ms, markerfacecolor=cbold[0], markeredgewidth=mew),
#                     # Line2D([0], [0], color=cfaded[1], lw=lw, label=labels[1], marker=marker, ms=ms, markerfacecolor=cbold[1], markeredgewidth=mew),
#                     # Line2D([0], [0], color=cfaded[2], lw=lw, label=labels[2], marker=marker, ms=ms, markerfacecolor=cbold[2], markeredgewidth=mew),
#                     Line2D([0], [0], color=cfaded[3], lw=lw, label=labels[0], marker=marker, ms=ms, markerfacecolor=cbold[3], markeredgewidth=mew),
#                     ]
#     fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(18.75,5.3), layout="constrained")


#     ###################
#     #    COLOR EVO    #
#     ###################
#     # This plots the whole timeseries
#     ax[0].plot(lightcurve[:,time_idx]/86400, color,
#                color=cfaded[3], lw=lw, marker=marker, ms=ms, fillstyle='none')

#     # This plots a single, filled point for the current snapshot
#     ax[0].plot(lightcurve[snapshot,time_idx]/86400, color[snapshot],
#                color=cbold[3], marker=marker, ms=ms, markeredgewidth=mew)

#     ax[0].set(title=r"Color",
#               xlabel="Days",
#               ylabel=r"$F_g$"+f"({g_center[-1]:.0f}) /" + r"$F_r$"f"({r_center[-1]:.0f}),
#               yscale="log")
#     # ax[0].set_yscale("log")
#     ax[0].legend(handles=legend_elements, title=f"time = {lightcurve[snapshot,time_idx]/86400:.1f} days", title_fontsize=14, fontsize=14)




    # plt.show()
    # break

    # figname = os.path.join(savedir, f"color_{snapshot:04}.png")
    # print("Saving figure:", figname)
    # plt.savefig(figname, format="png")
    # plt.close()
    # break
    
    # sys.exit()



