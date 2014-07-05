# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 221 - Adiabatenkoeffizient
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Autoren: Nils Fischer und Christian Kohlstedde
'''

import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import uncertainties as unc
import uncertainties.unumpy as unp

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats


# Konstanten
k_L_erw = 0.78 * 1.401 + 0.21 * 1.398 + 0.01 * 1.648
k_Ar_erw = 1.648
print k_L_erw, k_Ar_erw

#####
print u"# Methode nach Clement und Desormes"
#####

data = np.loadtxt('2.1.txt', skiprows=1)
d_h = 0.2
h1 = unp.uarray(data[:,0], d_h) - unp.uarray(data[:,1], d_h)
h3 = unp.uarray(data[:,2], d_h) - unp.uarray(data[:,3], d_h)
k = h1 / (h1 - h3)

print papstats.table(labels=['k'], columns=[k])

k = np.mean(k)
print papstats.pformat(k, label='k')
papstats.print_rdiff(k, k_L_erw)


#####
print u"# Methode nach Rüchhardt"
#####

p_0 = unc.ufloat(1004.6, 0.2) * const.hecto

def k_ruechhardt(m, V, r, n, t):
	global p_0
	T = t / n
	print papstats.pformat(T, 'T')
	p = p_0 + m * const.g / const.pi / r**2
	k = 4 * m * V / r**4 / T**2 / p
	return k

m = unc.ufloat(26.116, 0.002) * const.gram
V = unc.ufloat(5370, 5) * const.centi**3
r = 0.5 * unc.ufloat(15.95, 0.05) * const.milli
n = unc.ufloat(100, 1) 
t = unc.ufloat(60 + 37.26, 0.5)

k_L = k_ruechhardt(m, V, r, n, t)
print papstats.pformat(k_L, label='k_L')
papstats.print_rdiff(k_L, k_L_erw)

m = unc.ufloat(26.006, 0.002) * const.gram
V = unc.ufloat(5460, 5) * const.centi**3
r = 0.5 * unc.ufloat(15.97, 0.05) * const.milli
n = unc.ufloat(100, 1)
t = unc.ufloat(60 + 33.42, 0.5)

k_Ar = k_ruechhardt(m, V, r, n, t)
print papstats.pformat(k_Ar, label='k_Ar')
papstats.print_rdiff(k_Ar, k_Ar_erw)
