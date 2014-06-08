# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 233 - Fourieroptik
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
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


# Eichung

d = unp.uarray([2. * 71.3, 2. * 52.8, 2. * 31.0], 2 * 0.2) * const.milli / 100.
d_px = unp.uarray([1380.44, 1297.31, 1188.49], 5) - unp.uarray([710.88, 803.07, 908.87], 5)
px = np.mean(d / d_px)
print "Eichung: " + papstats.pformat(px / const.milli, unit='mm/px') + " <=> " + papstats.pformat(1. / px * const.milli, unit='px/mm')


#####
print "# 3.1: Quantitative Beobachtungen am Einzelspalt"
#####

# Messung bei hoher Intensität
n, x, dx, I, dI = np.loadtxt('2.1.2.txt', skiprows=1, unpack=True)
x = unp.uarray(x, dx)
I = unp.uarray(I, dI)
# Messung bei niedriger Intensität
n_low, x_low, dx_low, I_low, dI_low = np.loadtxt('2.1.1.txt', skiprows=1, unpack=True)
x_low = unp.uarray(x_low, dx_low)
I_low = unp.uarray(I_low, dI_low)

# Abstände vom 0. Maximum
i_0 = (len(n) - 1) / 2
i_0_low = (len(n_low) - 1) / 2
x_0 = x_low[i_0_low]
x = x - x_0
x[i_0] = unc.ufloat(0, 0)
x_low = x_low - x_0

# Untergrund
I = I - np.mean(unp.uarray([223.52, 220.72], 5))
I_low = I_low - np.mean(unp.uarray([156.4, 128.21], 5))

# Normierung
I_low = I_low / I_low[i_0_low]
I = I / np.mean([I[i_0+2], I[i_0-2]]) * np.mean([I_low[i_0_low+2], I_low[i_0_low-2]]) / I_low[i_0_low]
I[i_0] = I_low[i_0_low]

# Merge into combined data
for i_low in range(len(I_low)):
	i = i_0 - i_0_low + i_low
	x[i] = np.mean([x[i], x_low[i_low]])
	I[i] = np.mean([I[i], I_low[i_low]])

# Ordnung
n = (np.arange(len(I)) - i_0) / 2.
for i in range(len(n)):
	if i < i_0:
		n[i] = n[i] - 0.5
	elif i > i_0:
		n[i] = n[i] + 0.5

print papstats.table(labels=['n', 'x', 'I'], units=[None, 'px', None], columns=[n, x, I])

def fit_linear(x, m):
	return m * x

def approximate_I(n):
	x = n*const.pi + 1e-6
	return np.sin(x)**2/x**2


# Fit Minima
popt, pstats = papstats.curve_fit(fit_linear, n[1::2], x[1::2])

# Plot Abstand
plt.clf()
plt.title(u'Diagramm 3.1: Abstand der Interferenzmaxima und -minima vom Hauptmaximum')
papstats.plot_data(n[1::2], x[1::2], c='b', label='Minima')
papstats.plot_data(n[0::2], x[0::2], c='r', label='Maxima')
papstats.plot_fit(fit_linear, popt, xspace=n)
plt.legend(loc='upper left')
papstats.savefig_a4('3.1.png')

print papstats.pformat(80*const.milli * 635*const.nano / (popt[0] * px) / const.micro, label='Spaltweite d', unit='um')

# Maxima nach Fit
n_maxima_fit = x[0::2] / popt[0]
diff = np.abs(n[0::2]-n_maxima_fit)
print papstats.table(labels=['n_erw', 'n', u'|n-n_erw|', u'σ-Bereich'], columns=[n[0::2], n_maxima_fit, diff, unp.nominal_values(diff)/unp.std_devs(diff)])

# Plot Intensität
plt.clf()
plt.title(u'Diagramm 3.2: Normierte Intensitätsmaxima und -minima')
papstats.plot_data(n, I)
nspace = np.linspace(n[0], n[-1], 300)
Ispace = approximate_I(nspace)
plt.plot(nspace, Ispace, c='black')
papstats.savefig_a4('3.2.png')

# Vergleich der Intensitäten mit Erwartungswert
I_max = I[0::2]
I_max_erw = approximate_I(n[0::2])
diff = np.abs(I_max-I_max_erw)
print papstats.table(labels=['Ordnung', 'I', 'I_erw', u'|I-I_erw|', u'σ-Bereich'], columns=[np.arange(len(I_max))-(len(I_max)-1)/2., I_max, I_max_erw, diff, unp.nominal_values(diff)/unp.std_devs(diff)])







