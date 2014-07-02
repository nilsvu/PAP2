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
t = unp.uarray(np.mean(t, axis=1), np.std(t, axis=1)/np.sqrt(len(t)))

v = s / t
v_k = v / drho
y = (1 + 2.1 * r / R)
v_kl = v_k * y
r_sq = r**2

def fit_linear_origin(x, m):
	return m * x

popt, pstats = papstats.curve_fit(fit_linear_origin, r_sq, v_k)
popt_l, pstats_l = papstats.curve_fit(fit_linear_origin, r_sq, v_kl)

eta = 2. / 9. * const.g / popt_l[0]
print papstats.pformat(eta, label='eta')
v_lam = 2. / 9. * const.g * drho / eta * r_sq

plt.clf()
plt.title(u'Diagramm 3.1: Bestimmung der Viskosität nach Stokes')
papstats.plot_data(r_sq / const.centi**2, v_k, label='Messwerte')
papstats.plot_data(r_sq / const.centi**2, v_kl, label='Ladenburgkorrigierte Messwerte', color='red')
papstats.plot_data(r_sq / const.centi**2, v_lam / drho, label='Erwartungswerte', color='orange')
papstats.plot_fit(fit_linear_origin, popt, pstats, xspace=unp.nominal_values(r_sq), xscale=1./const.centi**2, eq=r'\frac{v}{\rho_K-\rho_F}=m*r^2', punits=[r'\frac{m^2}{kg*s}'])
papstats.plot_fit(fit_linear_origin, popt_l, pstats_l, xspace=unp.nominal_values(r_sq), xscale=1./const.centi**2, eq=r'\frac{v}{\rho_K-\rho_F}=m*r^2', punits=[r'\frac{m^2}{kg*s}'])
plt.xlabel('$r^2 \, [cm^2]$ mit $r$: Kugelradius')
plt.ylabel(r'$\frac{v}{\rho_K-\rho_F}$ mit $v$: mittlere Sinkgeschwindigkeit')
plt.legend(loc='upper left')
papstats.savefig_a4('1.png')

V = v / v_lam
Re = rho_F * v * 2. * r / eta

plt.clf()
plt.title(u'Diagramm 3.2: Abschätzung der kritischen Reynoldszahl')
plt.xscale('log')
papstats.plot_data(Re, V, color='black')
plt.xlabel('$Re$')
plt.ylabel(r'$\frac{v}{v_{lam}}$')
papstats.savefig_a4('2.png')

print papstats.table(labels=['d=2r', 'r^2', 'v', 'v/(rho_K-rho_F)', u'λ', 'v_k', 'v_k/(rho_K-rho_F)', 'v_lam', 'v/v_lam', 'Re'], units=['mm', 'mm^2', 'mm/s', 'mm^4/(g*s)', None, 'mm/s', 'mm^4/(g*s)', 'mm/s', None, None], columns=[2 * r / const.milli, r_sq / const.milli**2, v / const.milli, v_k / const.milli**4 * const.gram, y, v * y / const.milli, v_kl / const.milli**4 * const.gram, v_lam / const.milli, V, Re])

#####
print u"#2: Bestimmung der Viskosität nach Hagen-Poiseuille"
#####

Vt = unc.ufloat(25, 0.2 * np.sqrt(2)) * const.centi**3 / (9 * 60 + 39.08)
L = unc.ufloat(100, 0.5) * const.milli
R = unc.ufloat(1.5, 0.01) / 2. * const.milli
h_A, h_E = 536, 530.5

dp = rho_F * const.g * unc.ufloat(h_A + h_E, h_A - h_E) / 2. * const.milli
print dp, R, Vt, L

eta_2 = const.pi * dp * R**4 / 8. / Vt / L

papstats.print_rdiff(eta_2, eta)

print papstats.pformat(rho_F * Vt / const.pi / R * 2 / eta_2, label='Re')

