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

data = np.loadtxt('2.1.txt', skiprows=1)

d_h = [0.2 for i in range(len(data[:,0]))]

h1 = unp.uarray(data[:,0], d_h) - unp.uarray(data[:,1], d_h)

h3 = unp.uarray(data[:,2], d_h) - unp.uarray(data[:,3], d_h)

k = h1 / (h1 - h3)


p_0 = unc.ufloat(1004.6, 0.2) * const.hecto
Temperatur = unc.ufloat(23.3, 0.1) + const.zero_Celsius

def adiabtenkoeffizient_ruechhardt(m, V, r, n, t, p_0):
	T = t / n
	p = p_0 + m*const.g / const.pi / r**2
	k = 4 * m * V / r**4 / T**2 / p
	return k

m = unc.ufloat(26.116, 0.002) * const.gram
V = unc.ufloat(5370, 5) * const.centi**3
r = 0.5 * unc.ufloat(15.95, 0.05) * const.milli
n = unc.ufloat(100, 1) 
t = unc.ufloat(1*60 + 37.26, 0.5)

print adiabtenkoeffizient_ruechhardt(m, V, r, n, t, p_0)

m = unc.ufloat(26.006, 0.002) * const.gram
V = unc.ufloat(5460, 5) * const.centi**3
r = 0.5 * unc.ufloat(15.97, 0.05) * const.milli
n = unc.ufloat(100, 1)
t = unc.ufloat(1*60 + 33.42, 0.5)

print adiabtenkoeffizient_ruechhardt(m, V, r, n, t, p_0)