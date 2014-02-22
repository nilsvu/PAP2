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
papstats.plot_data(U, N)
papstats.savefig_a4('3.1.png')

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
print N1,N2
papstats.print_rdiff(N1[0],N2[0])
papstats.print_rdiff(N1[1],N2[1])

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

popt_gauss, pstats_gauss = papstats.curve_fit(fit_gauss, N[sl], n[sl], p0=[57.768, 7.659, 2094])
print 'Gauss:'
print popt_gauss
print pstats_gauss

# Poissonverteilung

def fit_poisson(x, m, A):
    return A*np.exp(-m)*(m**x)/gamma(x+1)

popt_poisson, pstats_poisson = papstats.curve_fit(fit_poisson, N[sl], n[sl], p0=[57.768, 2094])
print 'Poisson:'
print popt_poisson
print pstats_poisson

# plot
plt.clf()
papstats.plot_data(N, n)
papstats.plot_fit(fit_gauss, popt_gauss, pstats_gauss, np.linspace(0,100,num=200))
papstats.plot_fit(fit_poisson, popt_poisson, pstats_poisson, np.linspace(0,100,num=200))
xrange = 40
plt.xlim(int(popt_gauss[0].n-xrange),int(popt_gauss[0].n+xrange))
papstats.savefig_a4('3.2.png')

plt.clf()
plt.hist(fit_gauss(unp.nominal_values(N), *unp.nominal_values(popt_gauss))-n, bins=30)
plt.hist(pstats_gauss.residue, bins=30)
plt.hist(pstats_poisson.residue, bins=30)
papstats.savefig_a4('3.2.b.png')


#####
print('\n# 5 Vergleich der Poisson- und Gauß- Verteilung bei sehr kleinen Zählraten')
#####

t = 0.1

data = np.loadtxt('5.dat')

N = data[:,0]
N = unp.uarray(N,np.sqrt(N))
n = data[:,1]

sl = n > 10 # Häufigkeit mindestens 10

# Gaussverteilung

#popt_gauss, pstats_gauss = papstats.curve_fit(fit_gauss, N[sl], n[sl], p0=[4.134, 2.050, 5246])
print 'Gauss:'
print popt_gauss
print pstats_gauss

# Poissonverteilung

#popt_poisson, pstats_poisson = papstats.curve_fit(fit_poisson, N[sl], n[sl], p0=[4.134, 5246])
print 'Poisson:'
print popt_poisson
print pstats_poisson

# plot
plt.clf()
papstats.plot_data(N, n)
#papstats.plot_fit(fit_gauss, popt_gauss, pstats_gauss, np.linspace(0,100,num=200))
#papstats.plot_fit(fit_poisson, popt_poisson, pstats_poisson, np.linspace(0,100,num=200))
xrange = 40
#plt.xlim(int(popt_gauss[0].n-xrange),int(popt_gauss[0].n+xrange))
papstats.savefig_a4('3.3.png')


