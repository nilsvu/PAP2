# -*- coding: utf-8 -*-
__author__ = 'christian'

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c
import scipy.integrate
import scipy.optimize as opt
import sys
plt.ion()
plt.cla()

fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
fig.set_dpi(80)

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

# Lade Frequenzgang Daten
data = np.loadtxt('3.2.txt', skiprows=1)
f = data[:,0]
U_aus = data[:,1]


U_ein = 0.20
D = 0.001
dD_rel = 2e-2

g_f = 1 / D * U_aus / U_ein

def fit_g(f, V, o1, o2, n1, n2):
    return V / np.sqrt((1 + 1 / (f/o1)**(2*n1)) * (1 + (f/o2)**(2*n2)))


lower_f = 4.5e2
upper_f = 1e5
selectedRange = (lower_f < f) & (f < upper_f)
popt, pcov = opt.curve_fit(fit_g, f[selectedRange], g_f[selectedRange], p0=[1000, 1000, 50000, 5, 5])

def integrate_g_f(f, popt):
    return fit_g(f, *popt)**2

integrated = scipy.integrate.quad(integrate_g_f, 0, np.infty, args=popt)
B = integrated[0]
dB = integrated[1]

pstat = papstats.PAPStats(g_f[selectedRange], fit_g(f[selectedRange], *popt), ddof=5)
plt.xscale('log')
plt.yscale('log')
plt.plot(f, g_f,'.', ls='none', label='Messdaten')
plt.vlines(lower_f, 5, 1e3)
plt.vlines(upper_f, 5, 1e3)

plt.plot(f[selectedRange], fit_g(f[selectedRange], *popt), label=r'Fit von $g(f)$')
plt.title('Frequenzgang des Messaufbaus')
plt.xlabel(r'Frequenz $f$ in $Hz$')
plt.ylabel(r'$g(f)=\frac{1}{D} \frac{U_{aus}}{U_{ein}}$')

plt.legend()
text = r"$g(f)=\frac{V}{\sqrt{1+(\frac{\Omega_1}{f})^{2*n_1}}\sqrt{1+(\frac{f}{\Omega_2})^{2*n_2}}}$" + "\n"
parameters = np.array([
    ['V', None],
    [r"\Omega_1", 'Hz'],
    [r"\Omega_2", 'Hz'],
    [r"n_1", None],
    [r"n_2", None]
])

for i in range(parameters.shape[0]):
    text += papstats.pformat(v=popt[i], dv=pcov[i,i], label=parameters[i,0], unit=parameters[i,1]) + "\n"

text += pstat.legendstring().replace(', ', "\n") + "\n"
text += r"$B = \int_0^\infty g(f)^2 df$" + "\n"
text += papstats.pformat(v=B, dv=None, label='B', prec=3)

plt.text(2e3, 6e2, text, verticalalignment='top', backgroundcolor='#eeeeee')
plt.savefig('frequenzgang.png')

plt.cla()



data = np.loadtxt('2.1.txt', skiprows=1)

N = data[1:, 0]
R = data[1:, 1] * c.kilo
dR = R * 0.5e-2

U_aus = data[1:, 2] * c.milli #- 0.014073
dU_aus = np.max([data[1:, 3] * c.milli / np.sqrt(N), U_aus*3e-3], axis=0)

U_V = data[0, 2] * c.milli
dU_V = np.max([data[0, 3] / data[0,0], U_V*3e-2], axis=0)

U2 = U_aus**2 - U_V **2
dU2 = np.sqrt((2*U_aus*dU_aus)**2 + (2*U_V*dU_V)**2)

plt.errorbar(R, U2, xerr=dR, yerr=dU2, ls='none')

def fit_gerade(x, c):
    return  x*c

popt, pcov = opt.curve_fit(fit_gerade, R, U2, sigma=dU2)

x = np.linspace(0, 35e3)
plt.plot(x, fit_gerade(x, *popt))

T = 23.0 + c.zero_Celsius
k = popt[0]/ 4.0 / T / B