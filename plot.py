
import numpy as np
import pylab as plt

f = np.loadtxt(open('./data/intensity.dat','r'))
g = np.loadtxt(open('./data/grid.dat','r'))

U  = f[:,2]
V  = f[:,3]
Ip = f[:,4]
Im = f[:,5]
z  = g[:,1]

figure, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2,figsize=(12,6),sharex=True)

ax1.plot(z,U)
ax1.set_yscale('log')
ax1.set_title("Mean Intensity")
ax1.set_ylabel("U")

ax2.plot(z,V)
ax2.set_title("Flux")
ax2.set_ylabel("V")

ax3.plot(z,Ip)
ax3.set_yscale('log')
ax3.set_title("Upward Intensity")
ax3.set_xlabel("z")
ax3.set_ylabel("I+")

ax4.plot(z,Im)
ax4.set_yscale('log')
ax4.set_title("Downward Intensity")
ax4.set_xlabel("z")
ax4.set_ylabel("I-")

plt.tight_layout()
plt.show()
