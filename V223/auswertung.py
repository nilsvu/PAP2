# -*- coding: utf-8 -*-

'''
Auswertung des Versuches 233 - Brownsches Bewegung
Phyiskalisches Anfänger Praktikum II - Universität Heidelberg

Autoren: Christian Kohlstedde und Nils Fischer
'''

import numpy as np
import scipy as sp
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.constants as const
import scipy.interpolate as ip
import prettytable as pt

plt.ion()
plt.cla()

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

# Konstante Daten
T  = unc.ufloat(23.4 + sp.constants.zero_Celsius, 0.1)
a  = unc.ufloat(742, 15) * sp.constants.nano
eta = unc.ufloat(9.25e-4, 0.1e-4)

# Import der Daten
n, t, x, y = np.loadtxt('Messung.dat', skiprows=1, unpack=1)

plt.axis('equal')
plt.plot(x, y, marker='.', ms=10)
plt.xlabel(ur'x in $\mu m$')
plt.ylabel(ur'y in $\mu m$')
plt.title(ur'Abb.1: Bewegung eines Partikels')

papstats.savefig_a4('1.png')
plt.cla()

plt.axis('auto')

# Einheiten fürs Rechnen
x = x * sp.constants.micro
y = y * sp.constants.micro

# Erste Berechnungen
dx = np.diff(x)
dy = np.diff(y)
dt = np.diff(t)

dx2 = dx**2
dy2 = dy**2

r2 = dx2 + dy2


# Making the statistics table
labels = [u'Größe', 'Mittelwert', 'Standardabweichung', 'SE des Mittelwerts', 'Minimum', 'Maximum']

col1 = ['dt', 'dx', 'dy', 'dx^2', 'dy^2', 'r^2']

vars = np.array([dt, dx, dy, dx2, dy2, r2])

mean     = vars.mean(axis=1)
std      = vars.std(axis=1)
std_mean = std / np.sqrt(n.max())
min      = vars.min(axis=1)
max      = vars.max(axis=1)

columns = [col1, mean, std, std_mean, min, max]

table = pt.PrettyTable()

for i in range(len(labels)):
    table.add_column(labels[i], columns[i])

print table


# Werte aus der Tabelle
single_t = unc.ufloat(mean[0], std_mean[0])
single_r2 = unc.ufloat(mean[-1], std_mean[-1])


k = 6 * np.pi * eta * a * single_r2 / 4 / T / single_t
D = k * T / (6 * np.pi * eta * a)

print papstats.pformat(k, format='P', label='k', unit='J/K')
print papstats.pformat(D, format='P', label='D', unit='m^2/s')

# Histogram (Kombiniere und flachen des/der Arrays)
hist_data = (np.array([dx, dy]) / sp.constants.micro).flatten()

# Mathematische Werte
ma_mu, ma_sigma = hist_data.mean(), hist_data.std()

bin_count = 13
N, bins, patches = plt.hist(hist_data, bins=bin_count)
#N, bins = np.histogram(hist_data, bins=bin_count)

bin_middle = np.resize(bins, bin_count) + np.diff(bins)/2.0

# GaussFit
def fit_gauss(x, mu, sigma, A):
    return A/np.sqrt((sigma**2)*2*const.pi)*np.exp(-((x-mu)**2)/2/(sigma**2))

popt, pstats = papstats.curve_fit(fit_gauss, bin_middle, N)

mu, sigma = popt[0], popt[1]

xspace = np.linspace(-4, 4, 3000)
#y = mlab.normpdf(x, mu, sigma)
#plt.plot(x, y)

papstats.plot_fit(fit_gauss, popt, xspace=xspace)

plt.legend()

D = (sigma * sp.constants.micro)**2 / 2 / single_t
k = 6 * np.pi * eta * D * a / T

print papstats.pformat(k, format='P', label='k', unit='J/K')
print papstats.pformat(D, format='P', label='D', unit='m^2/s')

plt.xlabel(ur'$\Delta x$ und $\Delta y$ - Verschiebung in $\mu m$')
plt.ylabel(ur'$N$')
plt.title(ur'Abb.2: Histogramm der Verschiebungen mit gefitteter Gaußkurve')

papstats.savefig_a4('2.png')
plt.cla()

r2_cum = np.cumsum(r2)
new_t = np.resize(t, r2_cum.size) - 1

plt.plot(new_t, r2_cum)

def linear(x, m):
    return x * m

popt, pstats = papstats.curve_fit(linear, new_t, r2_cum)

xspace = np.linspace(0, 160, 2)
papstats.plot_fit(linear, popt, xspace=xspace)
plt.ylim((0, 3.5e-10))
plt.legend()

D = popt[0] / 4
k = 6 * np.pi * eta * D * a / T

print papstats.pformat(k, format='P', label='k', unit='J/K')
print papstats.pformat(D, format='P', label='D', unit='m^2/s')

plt.xlabel(ur'Zeit t in $s$')
plt.ylabel(ur'$r^2$ kumuliert')
plt.title(ur'Abb.3: Kumulative Verschiebung eines Partikels')

papstats.savefig_a4('3.png')