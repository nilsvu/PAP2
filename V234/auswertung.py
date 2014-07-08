# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 234 - Lichtquellen
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Autoren: Nils Fischer
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
print u"# 1: Sonnenspektrum"
#####

y, I_0 = np.loadtxt('1.1.txt', skiprows=2, unpack=True)
I_G = np.loadtxt('1.2.txt', skiprows=2)[:,1]

plt.clf()
plt.title(u'Diagramm 3.1: Sonnenspektrum bei geöffnetem und geschlossenem Fenster')
plt.plot(y, I_0, label=u'geöffnetes Fenster', color='black')
plt.plot(y, I_G, label=u'geschlossenes Fenster', color='grey')
plt.xlim(y[0], y[-1])
plt.xlabel(ur'Wellenlänge $\lambda \, [nm]$')
plt.ylabel(ur'Intensität $I \, [cnt]$')
plt.legend(loc='upper left')
papstats.savefig_a4('1.png')

# Absorption

A_G = 1 - I_G / I_0
y_min = 320
y_range = y > y_min

plt.clf()
plt.title(u'Diagramm 3.2: Absorption von Glas')
plt.plot(y[y_range], A_G[y_range], color='black')
plt.xlabel(ur'Wellenlänge $\lambda \, [nm]$')
plt.ylabel(ur'Absorption $A_G = 1 - \frac{I_G}{I_0}$')
plt.xlim(y_min, y[-1])
papstats.savefig_a4('2.png')

# Fraunhoferlinien

symbols, y_F, desc = np.loadtxt('fraunhoferlinien.txt', skiprows=1, delimiter='\t', dtype=object, unpack=True)
print papstats.table(labels=['Symbol', u'λ', u'Element / Molekül'], units=[None, 'nm', None], columns=[symbols, y_F, desc])
y_F = np.array(y_F, dtype=float)
dy_F = 1.5

y_range = (y > 350) & (y < 800)
I_min, I_max = 0, np.max(I_0) + 5000

plt.clf()
plt.title(u'Diagramm 3.3: Sonnenspektrum mit Fraunhoferlinien')
plt.plot(y[y_range], I_0[y_range], color='black')
axvl = None
for i in range(len(y_F)):
    axvl = plt.axvline(y_F[i], linestyle='dashed', lw=0.5, color='black')
    plt.fill_betweenx(np.linspace(I_min, I_max), y_F[i], y_F[i] + dy_F, color='grey', alpha=0.2, lw=0)
    ypos = 1000 + 4000 * i
    plt.annotate('$' + symbols[i] + '$', xy=(y_F[i], ypos), xytext=(y_F[i] + 10, ypos + 1500), arrowprops=dict(facecolor='black', width=0.5, frac=0.4, headwidth=4, shrink=0.1))
plt.ylim(I_min, I_max)
plt.xlabel(ur'Wellenlänge $\lambda \, [nm]$')
plt.ylabel(ur'Intensität $I \, [cnt]$')
fillbtwn = plt.Rectangle((0, 0), 1, 1, fc='grey', alpha=0.2, lw=0)
plt.legend([axvl, fillbtwn], [u'erwartete Fraunhoferlinien', ur'$\Delta \lambda = +' + str(dy_F) + u'nm$ geschätzte system. Abweichung'], loc='upper right')
papstats.savefig_a4('3.png')


#####
print "# 2: Natriumspektrum"
#####

y, I_low = np.loadtxt('3.1.txt', skiprows=1, unpack=True)
I_high = np.loadtxt('3.2.txt', skiprows=1)[:,1]

# Erwartete Linien

# Konstanten
E_Ry = -const.Rydberg * const.h * const.c / const.eV
hc = const.h * const.c / const.nano / const.eV

def fit_y_m(m, E_Ry, E_to, D):
    global hc
    return hc / (E_Ry / (m - D)**2 - E_to)

def fit_y_m_fixedERy(m, E_to, D):
    global hc, E_Ry
    return hc / (E_Ry / (m - D)**2 - E_to)

# 1. Nebenserie
y_m1 = unp.uarray([820.5, 569.5, 499.2, 467.9, 451.3, 440.1, 434.4, 0, 0, 0], [2, 1.5, 1.3, 1.1, 1.6, 1.5, 1.3, 0, 0, 0])
y_3p = y_m1[0]
E_3p = E_Ry / 3**2 - hc / y_3p
print papstats.pformat(E_3p, label='E_3p', unit='eV')
m_1 = np.arange(10) + 3
y_m1_erw = fit_y_m(m_1, E_Ry, E_3p, 0)

# 2. Nebenserie
y_m2 = unp.uarray([0, 617.1, 516.1, 476.3, 456.3], [0, 1.5, 1.3, 1.1, 1])
y_d = unc.ufloat(590, 4)
E_3s = E_3p - hc / y_d
print papstats.pformat(E_3s, label='E_3s', unit='eV')
d_s = 3 - unp.sqrt(E_Ry / E_3s)
print papstats.pformat(d_s, label='D_s')
m_2 = np.arange(5) + 4
y_m2_erw = fit_y_m(m_2, E_Ry, E_3p, d_s)

# 3. Nebenserie
y_m3 = unp.uarray([331.6, 0], [1.5, 0])
d_p = 3 - unp.sqrt(E_Ry / E_3p)
print papstats.pformat(d_p, label='D_p')
m_3 = np.arange(2) + 4
y_m3_erw = fit_y_m(m_3, E_Ry, E_3s, d_p)

def compare_y(m, y, y_erw):
    diff = y-y_erw
    diff[y==0] = 0
    print papstats.table(labels=['m', u'λ', u'λ_erw', u'Abweichung λ-λ_erw',  u'σ-Bereich'], units=[None, 'nm', 'nm', 'nm', None], columns=[m, y, y_erw, diff, np.abs(unp.nominal_values(diff)/unp.std_devs(diff))])

compare_y(m_1, y_m1, y_m1_erw)
compare_y(m_2, y_m2, y_m2_erw)
compare_y(m_3, y_m3, y_m3_erw)

# Plots

y_range = [(y > 300) & (y < 850), (y > 300) & (y < 540), (y > 550) & (y < 8500)]
I_offset = [0, 0, 0]

y_m_colors = ['blue', 'red', 'orange']
legend_loc = ['upper left', 'upper center', 'upper center']
title = [u'Natriumspektrum bei niedriger Intensität', u'Natriumspektrum bei hoher Intensität (Zoom)', u'Natriumspektrum bei hoher Intensität (Zoom)']

def plot_yvlines(y_m, y_m_erw, I_min, I_max, **kwargs):
    for i in range(len(y_m)):
        if y_m[i] != 0:
            plt.axvline(unp.nominal_values(y_m[i]), lw=0.5, **kwargs)
            plt.fill_betweenx(np.logspace(np.log10(I_min), np.log10(I_max)), y_m[i].nominal_value - y_m[i].std_dev, y_m[i].nominal_value + y_m[i].std_dev, alpha=0.2, lw=0, **kwargs)
        plt.axvline(unp.nominal_values(y_m_erw[i]), linestyle='dashed', lw=0.5, **kwargs)

for i in range(len(y_range)):
    plt.clf()
    plt.title(u'Diagramm 3.5.' + str(i+1) + ': ' + title[i])
    if i == 0:
        I = I_low
    else:
        I = I_high
    I = I[y_range[i]] + I_offset[i]
    I_min = np.min(I)
    I_max = np.max(I)
    plt.plot(y[y_range[i]], I, color='black')
    plot_yvlines(y_m1, y_m1_erw, I_min, I_max, color=y_m_colors[0])
    plot_yvlines(y_m2, y_m2_erw, I_min, I_max, color=y_m_colors[1])
    plot_yvlines(y_m3, y_m3_erw, I_min, I_max, color=y_m_colors[2])
    plt.axvline(y_d.nominal_value, lw=0.5, color='yellow')
    plt.fill_betweenx(np.logspace(np.log10(I_min), np.log10(I_max)), y_d.nominal_value - y_d.std_dev, y_d.nominal_value + y_d.std_dev, alpha=0.2, lw=0, color='yellow')
    plt.yscale('log')
    plt.xlabel(ur'Wellenlänge $\lambda \, [nm]$')
    plt.ylabel(ur'Intensität $I \, [cnt]$')
    plt.xlim(y[y_range[i]][0], y[y_range[i]][-1])
    plt.ylim(I_min, I_max)
    l = plt.Line2D((0, 0), (0, 0), color='grey', lw=0.5)
    l_fill = plt.Rectangle((0, 0), 1, 1, fc='grey', alpha=0.2, lw=0)
    l_erw = plt.Line2D((0, 0), (0, 0), color='grey', lw=0.5, linestyle='dashed')
    c1 = plt.Rectangle((0, 0), 1, 1, fc=y_m_colors[0], alpha=0.2, lw=0)
    c2 = plt.Rectangle((0, 0), 1, 1, fc=y_m_colors[1], alpha=0.2, lw=0)
    c3 = plt.Rectangle((0, 0), 1, 1, fc=y_m_colors[2], alpha=0.2, lw=0)
    plt.legend([l_erw, l, l_fill, c1, c2, c3], [u'erwartete Spektrallinie', u'Schätzwert der Line', u'Fehler der Schätzung', ur'1. Nebenserie: $md \rightarrow 3p$', ur'2. Nebenserie: $ms \rightarrow 3p$', ur'3. Nebenserie: $mp \rightarrow 3s$'], loc=legend_loc[i])
    papstats.savefig_a4('4.' + str(i+1) + '.png')


j = 1
def perform_fit(m, y_m, series):
    global j

    m = m[y_m!=0]
    y_m = y_m[y_m!=0]
    popt, pstats = papstats.curve_fit(fit_y_m, m, y_m, p0=[E_Ry, E_3p.nominal_value, 0])
    popt_fixedERy, pstats_fixedERy = papstats.curve_fit(fit_y_m_fixedERy, m, y_m, p0=[E_3p.nominal_value, 0])

    D_index = ['d', 's'][series-1]
    D = [0, d_s][series-1]

    plt.clf()
    plt.title(u'Diagramm 3.6.' + str(j) + u': Spektrallinien der ' + str(series) + u'. Nebenserie: $m' + D_index + ur' \rightarrow 3p$')
    papstats.plot_data(m, y_m)
    papstats.plot_fit(fit_y_m, popt, pstats, xspace=np.linspace(m[0], m[-1], 100), eq=r'\lambda_m = h*c/(E_{Ry}/(m - \Delta_' + D_index + ')^2 - E_{3p})', plabels=['E_{Ry}', 'E_{3p}', '\Delta_' + D_index], punits=['eV', 'eV', None])
    papstats.plot_fit(fit_y_m_fixedERy, popt_fixedERy, pstats_fixedERy, xspace=np.linspace(m[0], m[-1], 100), eq=r'E_{Ry}=13.605eV', plabels=['E_{3p}', '\Delta_' + D_index], punits=['eV', None])
    plt.ylabel(u'Wellenlänge der Spektrallinie $\lambda \, [nm]$')
    plt.xlabel(u'Anfangsniveau $m$')
    plt.legend(loc='upper right')
    papstats.savefig_a4('5.' + str(j) + '.png')
    j = j + 1

    erw = np.array([E_Ry, E_3p, D])
    res_fit = np.array(popt)
    diff = np.abs(erw - res_fit)
    res_fit_fixedERy = [E_Ry, popt_fixedERy[0], popt_fixedERy[1]]
    diff_fixedERy = np.abs(erw - res_fit_fixedERy)
    print papstats.table(labels=['', 'Erwartung', 'Fit', 'Abweichung', u'σ-Bereich', 'Fit mit E_Ry fixiert', 'Abweichung', u'σ-Bereich'], columns=[['E_Ry [eV]', 'E_3p [eV]', 'D_' + D_index], erw, res_fit, diff, unp.nominal_values(diff)/unp.std_devs(diff), res_fit_fixedERy, diff_fixedERy, unp.nominal_values(diff_fixedERy)/unp.std_devs(diff_fixedERy)])

perform_fit(m_1, y_m1, series=1)
perform_fit(m_2, y_m2, series=2)
