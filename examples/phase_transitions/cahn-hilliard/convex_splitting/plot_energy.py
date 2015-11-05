#!/usr/bin/python
import matplotlib.pylab as plt
import numpy as np

#dt = 0.00878906

t,iter,ener,mass = np.loadtxt('energy.log', unpack=True, skiprows=3, delimiter='\t')
#t = np.arange(0,len(ener))

#plt.plot(t, ener)
plt.semilogy(t, ener-np.min(ener))
plt.xlabel(r'Elapsed Time (taking exponential steps)')
#plt.xlabel(r'Time (using $\Delta t=0.009$, $Co=1.25$)')
plt.ylabel(r'$\mathcal{F}-\mathcal{F}_{min}=\sum\sum f_0\Delta x\Delta y$ (arb. units)')
plt.title("Periodic Cahn-Hilliard")
plt.savefig('energy.png', bbox_inches='tight', dpi=300)
plt.close()

plt.plot(t, mass)
#plt.semilogy(t, mass)
#plt.ylim([10, 5e4])
plt.xlabel(r'Elapsed Time (taking exponential steps)')
#plt.xlabel(r'Time (using $\Delta t=0.009$, $Co=1.25$)')
plt.ylabel(r'$m=\sum\sum c\Delta x\Delta y$ (arb. units)')
plt.title("Periodic Cahn-Hilliard")
plt.savefig('mass.png', bbox_inches='tight', dpi=300)
plt.close()


def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')
        
iter_av = movingaverage(iter, 75)

plt.plot(t, iter, 'k.')
plt.plot(t, iter_av, 'r')
plt.xlabel(r'Elapsed Time (taking exponential steps)')
#plt.xlabel(r'Time (using $\Delta t=0.009$, $Co=1.25$)')
plt.ylabel(r'Iterations to Converge (atol=$10^{-8}$)')
plt.title("Periodic Cahn-Hilliard")
plt.savefig('iter.png', bbox_inches='tight', dpi=300)
plt.close()

plt.plot(t, ener, label="Energy")
plt.plot(t, iter, label="Iterations")
plt.xlabel(r'Elapsed Time (taking exponential steps)')
#plt.xlabel(r'Time (using $\Delta t=0.009$, $Co=1.25$)')
plt.ylabel(r'Energy and Iterations')
plt.title("Periodic Cahn-Hilliard")
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
plt.savefig('compare.png', bbox_inches='tight', dpi=300)
plt.close()

print "%.1f iterations-per-timestep on average"%np.mean(iter)
