# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 212 - Zähigkeit von Flüssigkeiten
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


#####
print u"# 1: Bestimmung der Viskosität nach Stokes"
#####

t = np.loadtxt('1.txt')
r = np.array([9, 8, 7.144, 6, 5, 4, 3, 2, 1.5]) / 2. * const.milli * unc.ufloat(1, 0.01)
s = 30 * const.centi
T = unc.ufloat(25.25, 0.1) + 271.15
R = 75 / 2. * const.milli
rho_K_1 = np.array([1.36, 1.355, 1.375, 1.375, 1.375, 1.375, 1.375, 1.375, 1.375]) * const.gram/const.centi**3
rho_K_2 = np.array([1.365, 1.36, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38]) * const.gram/const.centi**3
rho_K = unp.uarray((rho_K_2 + rho_K_1) / 2., (rho_K_2 - rho_K_1) / 2.)
rho_F = unc.ufloat(1.1442, 0.0002) * const.gram/const.centi**3
drho = rho_K - rho_F

t = np.reshape(t, (len(r), 5))
t = unp.uarray(np.mean(t, axis=1), np.std(t, axis=1))

v = s / t
v_k = v / drho
v_kl = v_k * (1 + 2.1 * r / R)
r_sq = r**2

def fit_linear_origin(x, m):
	return m * x

popt, pstats = papstats.curve_fit(fit_linear_origin, r_sq, v_k)
popt_l, pstats_l = papstats.curve_fit(fit_linear_origin, r_sq, v_kl)

eta = 2. / 9. * const.g / popt_l[0]
print papstats.pformat(eta, label='eta')
v_lam = 2. / 9. * const.g * drho / eta * r_sq

plt.clf()
papstats.plot_data(r_sq / const.centi**2, v_k, label='Messwerte')
papstats.plot_data(r_sq / const.centi**2, v_kl, label='Ladenburgkorrigierte Messwerte', color='red')
papstats.plot_data(r_sq / const.centi**2, v_lam / drho, label='Erwartungswerte', color='orange')
papstats.plot_fit(fit_linear_origin, popt, pstats, xspace=unp.nominal_values(r_sq), xscale=1./const.centi**2, eq=r'\frac{v}{\rho_K-\rho_F}=m*r^2')
papstats.plot_fit(fit_linear_origin, popt_l, pstats_l, xspace=unp.nominal_values(r_sq), xscale=1./const.centi**2, eq=r'\frac{v}{\rho_K-\rho_F}=m*r^2')
plt.xlabel('$r^2 \, [cm^2]$ mit $r$: Kugelradius')
plt.ylabel(r'$\frac{v}{\rho_K-\rho_F}$ mit $v$: mittlere Sinkgeschwindigkeit')
plt.legend(loc='upper left')
papstats.savefig_a4('1.png')

print papstats.table(labels=['v', 'v_k', 'v_lam'], columns=[v, v_kl, v_lam])


#####
print u"#2: Bestimmung der Viskosität nach Hagen-Poiseuille"
#####

Vt = unc.ufloat(25, 0.2) * const.centi**3 / (9 * 60 + 39.08)
L = unc.ufloat(100, 0.5) * const.milli
R = unc.ufloat(1.5, 0.01) / 2. * const.milli

V = v / v_lam
Re = rho_F * v * 2. * r / eta

plt.clf()
plt.xscale('log')
papstats.plot_data(Re, V)
plt.xlabel('$Re$')
plt.ylabel(r'$\frac{v}{v_{lam}}$')
papstats.savefig_a4('2.png')

dp = rho_F * const.g * np.mean([unc.ufloat(536, 1), unc.ufloat(530.5, 1)]) * const.milli
print dp, R, Vt, L

eta_2 = const.pi * dp * R**4 / 8. / Vt / L

papstats.print_rdiff(eta_2, eta)

print papstats.pformat(rho_F * Vt / const.pi / R * 2 / eta_2, label='Re')

