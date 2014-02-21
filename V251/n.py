# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 251 - Statistik des radioaktiven Zerfalls
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as const
import scipy.optimize as opt
from scipy.special import gamma

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

#####
print('# 2 (Messung der Zählrohcharackteristik)')
#####

data = np.loadtxt('2.txt', skiprows=1)

U = unp.uarray(data[:,0], 10)
N = data[:,1]
N = unp.uarray(N, np.sqrt(N))

plt.clf()
plt.errorbar(unp.nominal_values(U), unp.nominal_values(N), xerr=unp.std_devs(U), yerr=unp.std_devs(N))
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.1.png', dpi=144)

#####
print('\n# 3 (Untersuchung des Plateauanstiegs)')
#####

data = np.loadtxt('3.txt', skiprows=1)

t = data[:,0]
N1 = data[:,1]
N1 = unp.uarray(N1,np.sqrt(N1))
N2 = data[:,2]
N2 = unp.uarray(N2,np.sqrt(N2))

Ndiff = N1-N2


#####
print('\n# 4 (Verifizierung der statistischen Natur des radioaktiven Zerfalls)')
#####

t = 0.5

data = np.loadtxt('4.dat')

N = data[:,0]
N = unp.uarray(N,np.sqrt(N))
n = data[:,1]

sl = n > 10 # Häufigkeit mindestens 10

# Gaussverteilung

def fit_gauss(x, m, s, A):
    return A/np.sqrt((s**2)*2*const.pi)*np.exp(-((x-m)**2)/2/(s**2))

popt, pcov = opt.curve_fit(fit_gauss, unp.nominal_values(N[sl]), unp.nominal_values(n[sl]), sigma=unp.std_devs(N[sl]), p0=[57.768, 7.659, 2094])
popt_gauss = unp.uarray(popt, np.sqrt(np.diagonal(pcov)))
pstats_gauss = papstats.PAPStats(unp.nominal_values(n[sl]), fit_gauss(unp.nominal_values(N[sl]), *unp.nominal_values(popt_gauss)), sigma=unp.std_devs(N[sl]), ddof=3)
print 'Gauss:'
print popt_gauss
print pstats_gauss

# Poissonverteilung

def fit_poisson(x, m, A):
    return A*np.exp(-m)*(m**x)/gamma(x+1)

popt, pcov = opt.curve_fit(fit_poisson, unp.nominal_values(N[sl]), unp.nominal_values(n[sl]), sigma=unp.std_devs(N[sl]), p0=[57.768, 2094])
popt_poisson = unp.uarray(popt, np.sqrt(np.diagonal(pcov)))
pstats_poisson = papstats.PAPStats(unp.nominal_values(n[sl]), fit_poisson(unp.nominal_values(N[sl]), *unp.nominal_values(popt_poisson)), sigma=unp.std_devs(N[sl]), ddof=2)
print 'Poisson:'
print popt_poisson
print pstats_poisson

# plot
plt.clf()
plt.errorbar(unp.nominal_values(N), unp.nominal_values(n), xerr=unp.std_devs(N), ls='none')
xspace = np.linspace(0, 100, num=200)
plt.plot(xspace, fit_gauss(xspace, *unp.nominal_values(popt_gauss)))
plt.plot(xspace, fit_poisson(xspace, *unp.nominal_values(popt_poisson)))
xrange = 40
plt.xlim(int(popt_gauss[0].n-xrange),int(popt_gauss[0].n+xrange))
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.2.png', dpi=144)




