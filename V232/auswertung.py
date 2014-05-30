# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 232 - Michelson-Interferometer
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


# Konstanten
y_erw = 546.07 * const.nano


#####
print u"3.1: Wellenlänge"
#####

dx = ( unp.uarray([11.2201, 11.426, 10.844], 2e-3 ) - unp.uarray([10.958, 11.169, 10.550], 2e-3 ) ) * const.milli / 5.
dm = unp.uarray([200, 201, 220], 5.)
y_list = 2 * dx / dm
y = np.mean(y_list)
print papstats.pformat(y / const.nano, label='y', unit='nm', format='c')
papstats.print_rdiff(y / const.nano, y_erw / const.nano)


#####
print u"3.2: Brechungsindex von Luft"
#####

a = unc.ufloat(50, 0.05) * const.milli
T = unc.ufloat(24.6 + const.zero_Celsius, 2)
T_0 = const.zero_Celsius
p_0 = const.atm / const.torr

p = unp.uarray([[-746, -592, -438, -284, -130], [-745, -595, -443, -287, -133], [-745, -597, -444, -290, -136]], 5)
dm = np.array([0, 10, 20, 30, 40])
m_0 = np.array([48, 49, 49])

def fit_dm(x, d):
    global m_0_i
    return d * x + m_0_i


plt.clf()
plt.suptitle(u'Diagramm 3.1: Vorbeigezogene Interferenzmaxima über Druck in der Küvette')
xspace = np.linspace(- const.atm / const.torr, 0)
ax = None
d_list = []
for i in range(len(p)):
    m_0_i = m_0[i]
    popt, pstats = papstats.curve_fit(fit_dm, p[i], dm)
    ax = plt.subplot(len(p), 1, i + 1, sharex=ax, sharey=ax)
    plt.title('Messung '+str(i+1))
    papstats.plot_data(p[i], dm, label='Messwerte')
    papstats.plot_fit(fit_dm, popt, xspace=xspace, eq='\Delta m=p*d+\Delta m_0', punits=['Torr^{-1}'])
    if i != len(p) - 1:
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        plt.xlabel(u'Druck $p \, [Torr]$')
    plt.ylabel(u'$\Delta m$')
    plt.xlim(xspace[0], xspace[-1])
    plt.ylim(0, 50)
    plt.legend(loc='upper left')
    d_list.append(popt[0])
papstats.savefig_a4('3.1.png')

print d_list
d = np.mean(d_list)
print papstats.pformat(d, label='d', unit='Torr^(-1)', format='c')
print papstats.pformat(d / const.torr / const.micro, label='d', unit='uPa^(-1)', format='c')

n_0_red = y_erw / 2. / a * d * p_0 * T / T_0
print papstats.pformat(n_0_red, label='n-1', format='c')
papstats.print_rdiff(n_0_red, 28e-5)


#####
print u"3.3: Kohärenzlängen"
#####

L = np.abs( unp.uarray([11.099, 10.878, 10.915, 10.937], [0.01, 0.005, 0.002, 0.002]) - unp.uarray([10.822, 11.042, 11.006, 10.988], [0.01, 0.005, 0.002, 0.002]) ) * const.milli / 5.
dy = np.array([5, 9, 16.2, 33]) * const.nano
L_erw = y_erw**2 / dy

def fit_linear(x, e):
    return e * x

popt, pstats = papstats.curve_fit(fit_linear, L_erw, L)

plt.clf()
plt.title(u'Diagramm 3.2: Vergleich der gemessenen Kohärenzlängen mit den Erwartungswerten')
xspace = np.linspace(L_erw[0], L_erw[-1])
papstats.plot_data(L_erw / const.micro , L / const.micro)
papstats.plot_fit(fit_linear, popt, xspace=xspace, xscale=1. / const.micro, yscale=1. / const.micro, eq='L=e*L_{erw}')
plt.legend(loc='upper left')
plt.xlabel(ur'Erwartungswert $L_{erw}=\frac{\lambda^2}{\Delta \lambda} \, [\mu m]$')
plt.ylabel(ur'Messwert $L=\frac{|s_0^\prime-s_0^{\prime \prime}|}{5} \, [\mu m]$')
papstats.savefig_a4('3.2.png')

print papstats.table(labels=[u'∆y [nm]', 'L [um]', 'L_erw [um]', u'|L-L_erw|', 'Sigmabereich'], columns=[dy/const.nano, L/const.micro, L_erw/const.micro, np.abs(L-L_erw), unp.nominal_values(np.abs(L-L_erw))/unp.std_devs(L-L_erw)])