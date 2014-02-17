# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c
import scipy.optimize as opt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

# Dämpfung
D = 0.001
dDrel = 0.2
dD = D * dDrel

#####
print "2 (Frequenzgang)"
#####

data = np.loadtxt('3.txt', skiprows=1)

Uein = 0.2
#dUein

f = data[:,0]
df = None
Uaus = data[:,1]
#dUaus

gf = Uaus/Uein/D
#dgf

# Fit Frequenzgang

def fit_gf(f, V, O1, O2, n1, n2):
    return V/np.sqrt(1+1/(f/O1)**(2*n1))/np.sqrt(1+(f/O2)**(2*n2))

popt, pcov = opt.curve_fit(fit_gf, f, gf, p0=[1000,1000,50000,5,5], sigma=df)

plt.clf()
plt.scatter(f, gf)
plt.plot(f, fit_gf(f, *popt))
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('2.png', dpi=144)
