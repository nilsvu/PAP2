# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 252 - Aktivierung von Indium und von Silber mit thermischen Neutronen
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as const

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

# Zählrohr
#r_Z =

#####
print "# 1 Nulleffekt"
#####

n0 = unc.ufloat(120,np.sqrt(120))/(5.*60.)


#####
print "\n# 2 Absorption von b-Strahlung in Aluminium"
#####

data = np.loadtxt('1.txt', skiprows=1)

x = data[:,0]*const.milli
N = data[:,2]
n = unp.uarray(N, np.sqrt(N))/data[:,1]
n0b = unc.ufloat(149,np.sqrt(149))/(5.*60.)

n = n-n0b

plt.clf()
plt.yscale('log')
papstats.plot_data(x, n)
papstats.savefig_a4('3.1.png')


#####
print "\n# 3 Absorption von y-Strahlung in Blei"
#####

data = np.loadtxt('2.txt', skiprows=1)

x = data[:,0]*const.milli
N = data[:,1]
n = unp.uarray(N, np.sqrt(N))/60

n = n-n0

def fit_damp(x, mu, n_0):
    return n_0*np.exp(-mu*x)

popt, pstats = papstats.curve_fit(fit_damp, x, n)

plt.clf()
papstats.plot_data(x, n)
papstats.plot_fit(fit_damp, popt, pstats, np.linspace(x[0],x[-1]))
plt.legend()
papstats.savefig_a4('3.2.png')


#####
print "\n# 5"
#####

data = np.loadtxt('4.txt', skiprows=1)

p1 = data[:,0]
p2 = data[:,1]
p = unp.uarray(p1+(p2-p1)/2.,(p2-p1)/2.)
N = data[:,2]
n = unp.uarray(N, np.sqrt(N))/60.

plt.clf()
papstats.plot_data(p, n)
papstats.savefig_a4('3.4.png')

