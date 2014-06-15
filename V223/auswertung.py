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
import scipy.constants as const
import prettytable as pt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

# Konstanten
T = unc.ufloat(23.4 + const.zero_Celsius, 0.1) # Raumtemperatur
a = unc.ufloat(742, 15) * const.nano / 2 # Partikelradius
eta = unc.ufloat(9.25e-4, 0.1e-4) # Viskosität
k_B_erw = 1.38065e-23

#####
print "3.1: Berechnung des mittleren Verschiebungsquadrats"
#####

n, t, x, y = np.loadtxt('Messung.dat', skiprows=1, unpack=True)
x = x * const.micro
y = y * const.micro

# Plot Partikelbewegung
plt.cla()
plt.axis('equal')
plt.plot(x / const.micro, y / const.micro, marker='o', ms=5, color='black')
plt.xlabel(ur'Koordinate $x \, [\mu m]$')
plt.ylabel(ur'Koordinate $y \, [\mu m]$')
plt.title(ur'Diagramm 3.1: Bewegung des beobachteten Latex-Partikels')
papstats.savefig_a4('1.png')
plt.axis('auto')

# Differenzen
dt = np.diff(t)
dx = np.diff(x)
dy = np.diff(y)
dx2 = dx**2
dy2 = dy**2
r2 = dx2 + dy2

# Spaltenstatistik
quantities = np.array([dt, dx / const.micro, dy / const.micro, dx2 / const.micro**2, dy2 / const.micro**2, r2 / const.micro**2])
mean = quantities.mean(axis=1)
std = quantities.std(axis=1)
std_mean = std / np.sqrt(n.max())
minimum = quantities.min(axis=1)
maximum = quantities.max(axis=1)
print papstats.table(labels=['dt', 'dx', 'dy', 'dx^2', 'dy^2', 'r^2'], units=['s', 'um', 'um', 'um^2', 'um^2', 'um^2'], columns=np.transpose([mean, std, std_mean, minimum, maximum]), rowlabels=['Mittelwert', 'Standardabw.', 'SE des MW', 'Minimum', 'Maximum'], prec=5)

# Mittelwerte
t_mean = unc.ufloat(mean[0], std_mean[0])
r2_mean = unc.ufloat(mean[-1], std_mean[-1]) * const.micro**2
print papstats.pformat(r2_mean / const.micro**2, label='r^2', unit='um^2', format='c')
print papstats.pformat(t_mean, label='t', unit='s', format='c')

# Boltzmannkonstante
k_B = 3./2. * const.pi * eta * a / T / t_mean * r2_mean
print papstats.pformat(k_B, format='c', label='k_B', unit='J/K')
papstats.print_rdiff(k_B, k_B_erw)

# Diffusionskoeffizient
D = k_B * T / (6 * const.pi * eta * a)
print papstats.pformat(D, format='c', label='D', unit='m^2/s')


#####
print "3.2: Kontrollverteilung"
#####

plt.cla()
plt.title(ur'Diagramm 3.2: Histogramm der Verschiebungen mit Gauß-Fit')

# kombiniertes Histogramm der x- und y-Daten
dr = np.append(dx, dy) / const.micro
n, bins, patches = plt.hist(dr, bins=13, label="Messungen")
bin_centers = bins[:-1] + np.diff(bins) / 2.

# Gauss Fit
def fit_gauss(x, mu, sigma, A):
    return A/np.sqrt((sigma**2)*2*const.pi)*np.exp(-((x-mu)**2)/2/(sigma**2))
popt, pstats = papstats.curve_fit(fit_gauss, bin_centers, n, p0=[dr.mean(), dr.std(), 1./np.sum(n)])

xspace = np.linspace(popt[0].nominal_value - 4 * popt[1].nominal_value, popt[0].nominal_value + 4 * popt[1].nominal_value, 100)
papstats.plot_fit(fit_gauss, popt, xspace=xspace, plabels=['\mu', '\sigma', 'A'], punits=['\mu m', '\mu m', None])
plt.xlim(xspace[0], xspace[-1])
plt.xlabel(ur'Verschiebung $\Delta x$ und $\Delta y$ $[\mu m]$')
plt.ylabel(ur'Häufigkeit $N$')
plt.legend()
papstats.savefig_a4('2.png')

# Berechnung der Konstanten
mu, sigma = popt[0] * const.micro, popt[1] * const.micro
D_fit = sigma**2 / 2. / t_mean
print papstats.pformat(D_fit, format='c', label='D_fit', unit='m^2/s')
papstats.print_rdiff(D_fit, D)
k_B_fit = 6 * const.pi * eta * D_fit * a / T
print papstats.pformat(k_B_fit, format='c', label='k_B_fit', unit='J/K')
papstats.print_rdiff(k_B_fit, k_B_erw)
papstats.print_rdiff(k_B_fit, k_B)


#####
print "3.3: Kumulative Verteilung der Verschiebungsquadrate"
#####

r2_cum = np.cumsum(r2)
t = t[:-1] - t[0]

def fit_lin(x, m):
    return x * m

popt, pstats = papstats.curve_fit(fit_lin, t, r2_cum)

plt.cla()
plt.title(ur'Diagramm 3.3: Kumulative Verschiebung des beobachteten Latex-Partikels')
plt.plot(t, r2_cum / const.micro**2)
xspace = np.linspace(0, t[-1])
papstats.plot_fit(fit_lin, popt, xspace=xspace, yscale=1/const.micro**2, punits=['m^2/s'], eq='m*t')
plt.xlim(xspace[0], xspace[-1])
plt.xlabel(ur'Zeit $t \, [s]$')
plt.ylabel(ur'$<r^2> \, [\mu m^2]$ kumulativ')
plt.legend(loc='upper left')
papstats.savefig_a4('3.png')

D_cum = popt[0] / 4
k_B_cum = 6 * const.pi * eta * D_cum * a / T
print papstats.pformat(k_B_cum, format='c', label='k', unit='J/K')
papstats.print_rdiff(k_B_cum, k_B_erw)
papstats.print_rdiff(k_B_cum, k_B)
print papstats.pformat(D_cum, format='c', label='D', unit='m^2/s')
papstats.print_rdiff(D_cum, D)
