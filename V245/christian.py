# -*- coding: utf-8 -*-
__author__ = 'christian'

import numpy as np
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as c
import scipy.integrate
import scipy.optimize as opt
import sys

plt.ion()
plt.cla()

import os

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)
import papstats


# Daten der Helmholzspule

r_H = (295 * c.milli) / 2
a_H = 147 * c.milli
n_H = 124

# Daten der Flachspule

n_F = 4e3
A_F = 41.7 * c.centi ** 2


# B-Feld der Helmholtzspule
def B_H(I):
    return c.mu_0 * 8 / np.sqrt(125) * I * n_H / r_H


#################################
print '# 1) Induktionsgesetz'

print u'a) in Abhängigkeit von f'
data = np.loadtxt('2.a.txt', skiprows=1)

f = unp.uarray(data[:, 0], data[:, 1])
omega = f * 2 * c.pi

U_ss = unp.uarray(data[:, 2], data[:, 3])
U_ind = U_ss / 2.

I = unc.ufloat(4, 0.05)

# U_ind(t) = -B*A*N*omega*sin(omega*t)
# U_ind_max = -B*A*N*omega
def U_ind_max(omega, B):
    return B * A_F * n_F * omega


U_ind_exp = U_ind_max(omega, B_H(I))

popt, pcov = opt.curve_fit(f=U_ind_max, xdata=unp.nominal_values(omega), ydata=unp.nominal_values(U_ind),
                           sigma=unp.std_devs(omega + U_ind))
B = unc.ufloat(popt[0], np.sqrt(pcov[0, 0]))
pstats = papstats.PAPStats(ydata=unp.nominal_values(U_ind),
                           ymodel=U_ind_max(unp.nominal_values(omega), B.nominal_value),
                           sigma=unp.std_devs(omega + U_ind), ddof=1)

plt.title(u'Induktionspannung in Abhängigkeit von der Rotationsfrequenz der Flachspule')
plt.xlabel(ur'Frequenz $f$ in $Hz$')
plt.ylabel(ur'Max. Induktionsspannung $U_{ind}$ in $V$')

plt.errorbar(x=unp.nominal_values(f), xerr=unp.std_devs(f), y=unp.nominal_values(U_ind), yerr=unp.std_devs(U_ind),
             ls='none', label='Messwerte')
label = "Fit-Werte mit \n"
label += ur"$Û_{ind} = B \ast A_F \ast n_F \ast \omega$" + "\n"
label += ur"$n_F=%d$, " % n_F
label += '$' + papstats.pformat(v=A_F, prec=3, label='A_F', unit='m^2') + ur'$, '
label += "$" + papstats.pformat(v=(B * 1000.), label='B', unit='mT') + "$\n"
label += pstats.legendstring()
plt.plot(unp.nominal_values(f), U_ind_max(unp.nominal_values(omega), B.nominal_value), label=label)
plt.errorbar(x=unp.nominal_values(f), xerr=unp.std_devs(f), y=unp.nominal_values(U_ind_exp),
             yerr=unp.std_devs(U_ind_exp), ls='none', label='Erwartungswerte')

plt.legend(loc=2, borderpad=1)

plt.savefig('2.a.png')
plt.cla()

print u'b) in Abhängigkeit vom Spulenstrom I'

data = np.loadtxt('2.b.txt', skiprows=1)

f = unc.ufloat(10.1, 0.1)
omega = f * 2 * c.pi

I = unp.uarray(data[:, 0], data[:, 1])

U_ss = unp.uarray(data[:, 2], data[:, 3])
U_ind = U_ss / 2.

U_ind_exp = B_H(I) * A_F * n_F * omega


def U_ind_max(I, c):
    return I * c


popt, pstats = papstats.curve_fit(U_ind_max, I, U_ind)

plt.xlabel('Spulenstrom $I$ in $A$')
papstats.plot_data(I, U_ind, label="Messpunkte")
papstats.plot_fit(U_ind_max, popt, pstats, unp.nominal_values(I), eq=u"Û_{ind} = c * I")
papstats.plot_data(I, U_ind_exp, label='Erwartungswerte')

plt.legend(borderpad=1)

papstats.savefig_a4('2.b.png')

print "Induktionspannung bei periodischem Feldstrom"

Omega = unc.ufloat(104, 1) * 2 * c.pi  # Kreisfrequenz der Wechselspannung

data = np.loadtxt('3.a.txt', skiprows=1)

a = unp.uarray(data[:, 0], 2) / 360. * 2 * c.pi
U_ind = unp.uarray(data[:, 1], data[:, 2]) / 2


def fit_cos(a, c, d):
    return c * np.abs(np.cos(a + d))


popt, pstats = papstats.curve_fit(fit_cos, a, U_ind)

plt.clf()
plt.subplot(121)
plt.xlabel(r'Winkel $\alpha$ in $\pi$')
plt.ylabel(ur'Induktionsspannung $Û_{ind}$ in $V$')
papstats.plot_data(a / c.pi, U_ind, label='Messpunkte')
papstats.plot_fit(fit_cos, popt, pstats, np.linspace(0, 2 * c.pi, num=100), xscale=1. / c.pi,
                  eq=ur'Û_{ind}=c*|cos(\alpha + \phi)|', plabels=['c', r'\phi'], punits=['V', 'rad'])
plt.xlim(0, 2)
plt.ylim(0, 1.6)
plt.legend(loc='upper center')
# Kreisplot
plt.subplot(122, polar=True)
papstats.plot_data(a, U_ind, capsize=0)
papstats.plot_fit(fit_cos, popt, pstats, np.linspace(0, 2 * c.pi, num=100))
plt.xlabel(ur'Winkel $\alpha$')
plt.ylabel(ur'Induktionsspannung $Û_{ind}$ in $V$')
fig = plt.gcf()
fig.suptitle(u'Induktionsspannung in Abhängigkeit vom Winkel der Flachspule')
papstats.savefig_a4('3.a.png')

plt.clf()

# b) induzierte und angelegte Spannung

a = unc.ufloat(0, 2)

data = np.loadtxt('3.b.txt', skiprows=1)

Omega = unp.uarray(data[:, 0], data[:, 1]) * 2 * c.pi
U_H = unp.uarray(data[:, 2], data[:, 3]) / 2.
I_H = unp.uarray(data[:, 4], data[:, 5]) * c.milli
U_ind = unp.uarray(data[:, 6], 0.05) / 2.

V = U_ind / U_H

plt.title(u'Verhältnis von induzierter und angelegter Spannung als Funktion der Frequenz')
papstats.plot_data(Omega / c.kilo, V, label='Messpunkte')
plt.xlabel(ur'Wechselstromfrequenz $\Omega$ in $kHz$')
plt.ylabel(ur'Verhältnis $\frac{Û_{ind}}{Û_H}$')
plt.legend()
papstats.savefig_a4('3.b.png')

# c) Widerstand

R = U_H / I_H


def fit_R(Omega, L):
    return L * Omega


popt, pstats = papstats.curve_fit(fit_R, Omega, R)

plt.clf()
plt.title(u'Widerstand der Helmholtzspulen als Funktion der Frequenz')
papstats.plot_data(Omega / c.kilo, R, label='Messpunkte')
papstats.plot_fit(fit_R, popt, pstats, np.linspace(Omega[0].n, Omega[-1].n), xscale=1. / c.kilo, eq=r'Z=L*\Omega',
                  punits=['H'])
plt.xlabel(ur'Wechselstromfrequenz $\Omega$ in $kHz$')
plt.ylabel(ur'Widerstand $Z=\frac{Û_H}{I_H}$ in $\Omega$')
plt.legend(loc='lower right')
papstats.savefig_a4('3.c.png')
