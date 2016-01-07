#!/usr/bin/python
import matplotlib.pylab as plt
import numpy as np

dt = 0.0253125

iter,ener,mass = np.loadtxt('energy.log', unpack=True, skiprows=3, delimiter='\t')
t = dt*np.arange(0,len(ener))

#plt.plot(t, ener)
plt.semilogy(t, ener-np.min(ener))
plt.xlabel(r'Time (using $\Delta t=0.025$, Co=920)')
plt.ylabel(r'$\mathcal{F}-\mathcal{F}_{min}=\sum\sum f_0\Delta x\Delta y$ (arb. units)')
plt.title("Periodic Cahn-Hilliard")
plt.savefig('energy.png', bbox_inches='tight', dpi=300)
plt.close()

plt.semilogy(t, mass)
plt.xlabel(r'Time (using $\Delta t=0.025$, Co=920)')
plt.ylabel(r'$m=\sum\sum c\Delta x\Delta y$ (arb. units)')
plt.title("Periodic Cahn-Hilliard")
plt.savefig('mass.png', bbox_inches='tight', dpi=300)
plt.close()

print "%d total, %.1f iterations-per-timestep on average"%(np.sum(iter),np.mean(iter))
