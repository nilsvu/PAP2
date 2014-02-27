# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 252 - Aktivierung von Indium und von Silber mit thermischen Neutronen
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import scipy as sp
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as const
import scipy.interpolate as ip
import datetime as dt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

# Zählrohr
r_Z = 7*const.milli


#####
print "# 1 Nulleffekt"
#####

n0 = unc.ufloat(120,np.sqrt(120))/(5.*60.)


#####
print "\n# 2 Absorption von b-Strahlung in Aluminium"
#####

data = np.loadtxt('1.txt', skiprows=1)

x = unp.uarray(data[:,0]*const.milli, 0.05*const.milli)
N = data[:,2]
n = unp.uarray(N, np.sqrt(N))/data[:,1]
n0b = unc.ufloat(149,np.sqrt(149))/(5.*60.)

n = n-n0b

sl = (n > 0)

plt.clf()
plt.yscale('log')
plt.xlabel('Absorberdicke $x \, [mm]$')
plt.ylabel(u'korrigierte Zählrate '+r'$n-n_0^{\beta}$')
papstats.plot_data(x/const.milli, n)
xspace = np.linspace(0,0.004,num=1000)
for k in range(1,6):
    extrapolater = ip.UnivariateSpline(unp.nominal_values(x[sl]), np.log(unp.nominal_values(n[sl])), w=1./np.log(unp.std_devs(x+n)[sl]), k=k)
    plt.plot(xspace/const.milli, np.exp(extrapolater(xspace)), label='$k='+str(k)+'$', alpha=0.3)
sl_upper = 10
sl = sl & (n <= sl_upper)
k = 5
extrapolater = ip.UnivariateSpline(unp.nominal_values(x[sl]), np.log(unp.nominal_values(n[sl])), k=k)
plt.plot(xspace/const.milli, np.exp(extrapolater(xspace)), label='$k='+str(k)+'$ mit $n \in (0,'+str(sl_upper)+']$')
extrapolater = ip.UnivariateSpline(unp.nominal_values(x[sl]), np.log(unp.nominal_values(n[sl])+unp.std_devs(n[sl])), k=k)
plt.plot(xspace/const.milli, np.exp(extrapolater(xspace)), label='Fehlerkurve zu $k='+str(k)+'$ mit $n \in (0,'+str(sl_upper)+']$')
plt.ylim(1e-2,1e2)
plt.legend(loc='lower left', title='Extrapolation mit Polynomen\nder Ordnung:')
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
ax1 = plt.subplot(211)
papstats.plot_data(x, n)
papstats.plot_fit(fit_damp, popt, pstats, np.linspace(x[0],x[-1]), plabels=[r'\mu','n_0'])
plt.legend()
plt.subplot(212, sharex=ax1)
plt.yscale('log')
papstats.plot_data(x, n)
papstats.plot_fit(fit_damp, popt, pstats, np.linspace(x[0],x[-1]))
plt.setp(ax1.get_xticklabels(), visible=False)
papstats.savefig_a4('3.2.png')


#####
print "\n# 4 Bestimmung der Aktivität"
#####

# Erwartungswert
A0 = unc.ufloat(3.7e6,0)
t = (dt.date(2014,2,27)-dt.date(2012,2,2)).days
Th = 5.27*365.242199
A_erw = A0*np.exp(-np.log(2)/Th*t)

# Messdaten
x, dx, N = np.loadtxt('3.txt', skiprows=1, unpack=True)
x = unp.uarray(x, dx)
x += unc.ufloat(60,0.5)-unc.ufloat(160,0.5)+unc.ufloat(119,0.5)
x *= const.milli
N = unp.uarray(N, np.sqrt(N))
n = N/2./60.

# Berechnung der Aktivität
e = 0.04
A = 4*n/e*x**2/r_Z**2

# Vergleich mit Erwartungswerten
for i in range(len(A)):
    papstats.print_rdiff(A[i], A_erw)

# Raumwinkelkorrektur
l = 0.04
k1 = (4*n/e*(x+l/4.)**2/r_Z**2)/A
print k1, A*k1

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

