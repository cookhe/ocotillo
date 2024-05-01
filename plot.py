
import numpy as np
import pylab as plt

f = np.loadtxt(open('./intensity.dat','r'))

U  = f[:,2]
V  = f[:,3]
Ip = f[:,4]
Im = f[:,5]

figure, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2,figsize=(12,6),sharex=True)

ax1.plot(U)
ax1.set_yscale('log')

ax2.plot(V)

ax3.plot(Ip)
ax3.set_yscale('log')

ax4.plot(Im)
ax4.set_yscale('log')

plt.show()
