# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 213 - Kreisel
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Authoren: Nils Fischer und Christian Kohlstedde
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
print u"# 3.2: Dämpfung des Kreisels"
#####

def fit_exp(x, a, A):
	return A * np.exp(a * x)

t, w = np.loadtxt('2.2.txt', skiprows=1, unpack=True)
t = t
w = unp.uarray(w, 10)
popt, pstats = papstats.curve_fit(fit_exp, t, w)

plt.clf()
plt.title(u'Diagramm 3.1: Dämpfung des Kreisels')
papstats.plot_data(t, w)
papstats.plot_fit(fit_exp, popt, pstats, xspace=np.linspace(t[0], t[-1], 100), eq=r'w_0*e^{\lambda * t}', plabels=[r'\lambda', '\omega_0'], punits=[ur'\frac{1}{min}', ur'\frac{2π}{min}'])
plt.xlabel('Zeit $t \, [min]$')
plt.ylabel(ur'Drehfrequenz $\omega_F \, [\frac{2π}{min}]$')
plt.yscale('log')
plt.legend()
papstats.savefig_a4('3.1.png')

tau = -1. / popt[0] * const.minute
t_H = np.log(2) * tau
print papstats.pformat(tau / const.minute, label=u'Dämpfungskonstante tau', unit='min')
print papstats.pformat(t_H / const.minute, label=u'Halbwertszeit t_H', unit='min')


#####
print u"# 3.3: Präzession"
#####

data = np.loadtxt('2.3.txt', skiprows=1)
w_1 = unp.uarray(data[:,0], 10) * 2 * const.pi / const.minute
T_P = unp.uarray(np.transpose(data[:,1::2]), 2)
w_2 = unp.uarray(np.transpose(data[:,2::2]), 10) * 2 * const.pi / const.minute
w_2_erw = w_1 * unp.exp(-T_P / tau)
diff = w_2 - w_2_erw
w_2 = np.mean([w_2, w_2_erw], axis=0)
w_F = unp.uarray(unp.nominal_values(w_1 + w_2) / 2., np.sqrt((unp.nominal_values(w_1 - w_2) / 2)**2 + unp.std_devs(w_1 + w_2)**2))

def fit_linear_origin(x, m):
	return m * x

n_m = np.array([1, 1, 2, 2])
d = unp.uarray([15, 20, 15, 20], 0.1) * const.centi
mlabels = [str(n_m[i]) + 'm, ' + str(int(d[i].nominal_value/const.centi)) + 'cm' for i in range(len(n_m))]
d = d + 1.1 / 2. * const.centi * n_m
s = []
plt.clf()
plt.suptitle(u'Diagramm 3.2: Präzessionszeit über der Drehfrequenz verschiedener Drehmomente')
ax = None
for i in range(len(T_P)):
	print "Messung mit Drehmoment " +  mlabels[i] + ":"
	print papstats.table(labels=['w_1', 'T_P', 'w_2', 'w_2_erw', 'w_2-w_2_erw', u'σ-Bereich', 'MW w_F'], units=[u'2π/min', 's', u'2π/min', u'2π/min', u'2π/min', None, 'Hz'], columns=[w_1/2./const.pi*const.minute, T_P[i], w_2[i]/2./const.pi*const.minute, w_2_erw[i]/2./const.pi*const.minute, diff[i]/2./const.pi*const.minute, np.abs(unp.nominal_values(diff[i]))/unp.std_devs(diff[i]), w_F[i]])
	popt, pstats = papstats.curve_fit(fit_linear_origin, w_F[i], T_P[i])
	s.append(popt[0])
	ax_i = plt.subplot(2, 2, i + 1, sharex=ax)
	if ax is None:
		ax = ax_i
	plt.title('Messung mit Drehmoment ' + mlabels[i])
	papstats.plot_data(w_F[i], T_P[i])
	papstats.plot_fit(fit_linear_origin, popt, xspace=unp.nominal_values(w_F[i]), plabels=['s'], punits=['s^2'])
	if i >= 2:
		plt.xlabel('mittlere Drehfrequenz $\omega_F \, [Hz]$')
	else:
		plt.setp(ax_i.get_xticklabels(), visible=False)
	if i%2 == 0:
		plt.ylabel(u'Präzessionszeit $T_P \, [s]$')
	plt.legend(loc='upper left')
papstats.savefig_a4('3.2.png')

s = np.array(s)

I = n_m * 9.85 * const.gram * const.g * d * s / 2 / const.pi
print papstats.table(labels=['Messung', 's', 'd', 'I'], columns=[mlabels, s, d / const.centi, I / const.gram ], units=[None, 's^2', 'cm', 'g*m^2'])
I = np.mean(I)
print papstats.pformat(I / const.gram, unit='g*m^2', label=u'mittleres Trägheitsmoment')


#####
print u"# 3.4: Umlauf der momentanen Drehachse um die Figurenachse"
#####

w_F, t = np.loadtxt('2.4.txt', skiprows=1, unpack=True)
w_F = unp.uarray(w_F, 5) * 2 * const.pi / const.minute
W = 2 * const.pi / unp.uarray(t, 1) * 10

print papstats.table(labels=['w_F', u'Ω'], units=['Hz', 'Hz'], columns=[w_F, W])

popt, pstats = papstats.curve_fit(fit_linear_origin, w_F, W)

plt.clf()
plt.title(u'Diagramm 3.3: Umlauffrequenz der momentanen Drehachse über der Eigendrehfrequenz')
papstats.plot_data(w_F, W)
papstats.plot_fit(fit_linear_origin, popt, xspace=unp.nominal_values(w_F), plabels='s')
plt.xlabel('Eigendrehfrequenz $\omega_F \, [Hz]$')
plt.ylabel('Umlauffrequenz der momentanen Drehachse $\Omega \, [Hz]$')
plt.legend(loc='upper left')
papstats.savefig_a4('3.3.png')

dI = I / (1 / popt[0] - 1)
print papstats.pformat(dI / const.gram, label='I_x-I_z', unit='g*m^2')
I_x = dI + I
print papstats.pformat(I_x / const.gram, label='I_x', unit='g*m^2')


#####
print u"# 3.5: Nutation"
#####

w_F, w_N = np.loadtxt('2.5.txt', skiprows=1, unpack=True)
w_F = unp.uarray(w_F, 5)
w_N = unp.uarray(w_N, 15)

popt, pstats = papstats.curve_fit(fit_linear_origin, w_F, w_N)

plt.clf()
plt.title(u'Diagramm 3.4: Nutationsfrequenz über der Eigendrehfrequenz')
papstats.plot_data(w_F, w_N)
papstats.plot_fit(fit_linear_origin, popt, xspace=unp.nominal_values(w_F), plabels='s')
plt.xlabel('Eigendrehfrequenz $\omega_F \, [Hz]$')
plt.ylabel('Nutationsfrequenz $\Omega \, [Hz]$')
plt.legend(loc='upper left')
papstats.savefig_a4('3.4.png')

I_x_2 = I / popt[0]
print papstats.pformat(I_x_2 / const.gram, label='I_x_2', unit='g*m^2')
papstats.print_rdiff(I_x_2, I_x)







