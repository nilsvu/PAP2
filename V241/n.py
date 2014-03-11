# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 241 - Wechselstromeigenschaften von RLC-Gliedern
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import scipy as sp
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as const
import scipy.interpolate as ip
import prettytable as pt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats


#####
print "# 1: Bestimmung der Zeitkonstante eines RC-Glieds"
#####

C, R, f, T_H = np.loadtxt('1.txt', skiprows=1, converters=dict.fromkeys([3], unc.ufloat_fromstr), dtype=object, unpack=True)
C, R, f, T_H = np.array(C, dtype=float)*const.nano, np.array(R, dtype=float)*const.kilo, np.array(f, dtype=float), T_H*const.micro
C, R = unp.uarray(C, C*0.1), unp.uarray(R, R*0.05)

tau_exp = T_H/np.log(2)
tau_theo = R*C


print papstats.table(labels=['C', 'R', 'f', u'τ_exp', u'τ_theo', 'Abweichung', 'rel. Abw.', u'σ-Bereich'], units=['nF', u'kΩ', 'Hz', u'µs', u'µs', u'µs', None, None], columns=[C/const.nano, R/const.kilo, f, tau_exp/const.micro, tau_theo/const.micro] + list(papstats.rdiff(tau_exp/const.micro, tau_theo/const.micro)))

print "Messung aus Stromverlauf:", papstats.pformat(unc.ufloat(37.2,2)/np.log(2), unit='us')


#####
print "# 3: Frequenz- und Phasengang eines RC-Glied"
#####

# Theoretische Grenzfrequenz

R = 1*const.kilo * unc.ufloat(1, 0.05)
C = 47*const.nano * unc.ufloat(1, 0.1)
f_G_theo = 1./2./const.pi/R/C

print "Theoretische Grenzfrequenz:", papstats.pformat(f_G_theo/const.kilo, unit='kHz', format='c')

# Phasengang

f, dt = np.loadtxt('3.txt', skiprows=1, converters=dict.fromkeys([1], unc.ufloat_fromstr), dtype=object, unpack=True)
f, dt = np.array(f, dtype=float)*const.kilo, dt*const.micro
phi = 360*dt*f

plt.clf()
plt.title(u'Diagramm 3.1: Phasengang eines RC-Glieds')
plt.xlabel('Wechselstromfrequenz der Eingangsspannung $f \, [kHz]$')
plt.ylabel(ur'Phasenverschiebung der Ausgangsspannung $\Phi \, [^\circ]$')
plt.xscale('log')
xspace = np.logspace(0.2,1)*const.kilo
yspace = np.linspace(0,90)
plt.ylim(yspace[0], yspace[-1])
plt.xlim(xspace[0]/const.kilo, xspace[-1]/const.kilo)
papstats.plot_data(f/const.kilo, phi, label='Messpunkte')
extrapolater = ip.UnivariateSpline(unp.nominal_values(f), unp.nominal_values(phi), k=3)
plt.plot(xspace/const.kilo, extrapolater(xspace), label='Extrapolation', color='b')
plt.axhline(45, color='black', label='$\Phi = 45^\circ$')
plt.axvline(f_G_theo.n/const.kilo, label='Theoretische Grenzfrequenz $f_G^{theo}$\nmit Fehlerbereich', color='g')
plt.fill_betweenx(yspace, (f_G_theo.n-f_G_theo.s)/const.kilo, (f_G_theo.n+f_G_theo.s)/const.kilo, color='g', alpha=0.2)
f_G_phas = unc.ufloat(3.09, 0.2)*const.kilo
plt.axvline(f_G_phas.n/const.kilo, color='r', label='Grenzfrequenz $'+papstats.pformat(f_G_phas/const.kilo, unit='kHz', label='f_G^{phas}')+'$\nnach Extrapolation mit Fehlerbereich')
plt.fill_betweenx(yspace, (f_G_phas.n-f_G_phas.s)/const.kilo, (f_G_phas.n+f_G_phas.s)/const.kilo, color='r', alpha=0.2)
plt.legend()
papstats.savefig_a4('3.1.png')


#####
print "# 4: Frequenzgang eines Serienschwingkreises"
#####

data = np.loadtxt('4.txt', skiprows=1, converters=dict.fromkeys(range(1,6), unc.ufloat_fromstr), dtype=object)
R, f_R, df, U_E, U_A = np.array(data[:,0], dtype=float), data[:,1]*const.kilo, (data[:,3]-data[:,2])*const.kilo, data[:,4], data[:,5]
R = R * unc.ufloat(1, 0.05)

C = 47*const.nano
L = 1./(f_R*2*const.pi)**2/C

# Gesamtwiderstand
R_G = df*2*const.pi*L

# Verlustwiderstand
R_V_band = R_G-R
R_V_res = R*(U_E/U_A-1)

print u"Induktivität der Spule:"
print papstats.table(labels=['R', 'f_R', 'L', u'∆f', 'R_G', 'R_V^band', 'R_V^res'], units=[u'Ω', 'kHz', 'mH', 'kHz', u'Ω', u'Ω', u'Ω'], columns=[R, f_R/const.kilo, L/const.milli, df/const.kilo, R_G, R_V_band, R_V_res])

L = np.mean(L)
L1_1 = L
print u"Induktivität der Spule (gemittelt):", papstats.pformat(L/const.milli, format='c')
print u"Verlustwiderstand (gemittelt):", papstats.pformat(np.mean(np.append(R_V_band, R_V_res)), format='c')


#####
print u"# 5: Bestimmung der Dämpfungskonstanten"
#####

C = 47*const.nano * unc.ufloat(1, 0.1)
A = np.loadtxt('5.1.txt', skiprows=1, converters={0:unc.ufloat_fromstr}, dtype=object)
D = np.array([unc.umath.log(A[i]/A[i+1]) for i in range(len(A)-1)])
T = np.loadtxt('5.2.txt', skiprows=1, converters={0:unc.ufloat_fromstr}, dtype=object)*const.milli

print papstats.table(labels=['A_1', 'A_2', u'Λ', 'T'], units=['V', 'V', None, u'µs'], columns=[A[:-1], A[1:], D, T/const.micro])

D = np.mean(D)
T = np.mean(T)
print "Schwingungsperiode (gemittelt):", papstats.pformat(T/const.micro, unit=u'µs', label='T', format='c')
print "Log. Dekrement (gemittelt):", papstats.pformat(D, label='D', format='c')

d = D/T
L = 1./((2*const.pi/T)**2+d**2)/C
L1_2 = L
R_G_dekr = d*2*L

print u"Dämpfungskonstante:", papstats.pformat(d, label='d', format='c')
print u"Induktivität der Spule:", papstats.pformat(L/const.milli, label='L', unit='mH', format='c')
print "Verlustwiderstand:", papstats.pformat(R_G_dekr, format='c')

L1 = np.mean([L1_1, L1_2])

#####
print u"# 6: Resonanzüberhöhung"
#####

C = 47*const.nano * unc.ufloat(1, 0.1)
R = 220 * unc.ufloat(1, 0.1)
L = L1

f_R_exp = unp.uarray([4.04, 3.87, 4.19], 0.06)*const.kilo
f_R_theo = np.array([1./(2*const.pi*unc.umath.sqrt(L*C)), unc.umath.sqrt(1./L/C-2*d**2)/2./const.pi, unc.umath.sqrt(1./L/C+2*d**2)/2./const.pi])

print papstats.table(labels=['f_R^exp', 'f_R^theo', 'Abweichung', 'rel. Abw.', u'σ-Bereich'], units=['kHz', 'kHz', 'kHz', None, None], columns=[f_R_exp/const.kilo, f_R_theo/const.kilo]+list(papstats.rdiff(f_R_exp/const.kilo, f_R_theo/const.kilo)))


#####
print u"# 7: Bandsperre"
#####

f_R_theo = 1./unc.umath.sqrt(L*C)/2./const.pi
f_R_exp = unc.ufloat(4.06, 0.11)*const.kilo

papstats.print_rdiff(f_R_exp/const.kilo, f_R_theo/const.kilo)


#####
print u"# 8: Signalformung"
#####

f = np.array([0.103, 3.6, 7.93])*const.kilo
L_U = np.array([[-2.81, -11.25, -19.69], [-31.88, -13.75, -20.63], [-2.81, -15.31, -28.75], [-2.65, 2.98, -29.21], [-31.4, -12.02, -23.27], [-60.31, -23.44, -47.5]])
labels = ['Ungefiltert', 'RC-Hochpass', 'RC-Tiefpass', 'LC-Tiefpass', u'RLC-Bandpass $R=1kΩ$', u'RLC-Bandpass $R=47Ω$']

U = 10**(L_U/20.)
V = U/U[0]

tablestrings = []
rowlabels = ['L_U [dBV]', 'U [V/V_rms]', 'U/U_0']
for i in range(len(labels)):
    table = pt.PrettyTable(float_format='.2')
    table.add_column(labels[i], rowlabels)
    for j in range(len(f)):
        table.add_column('f_'+str(j+1), [L_U[i,j], U[i,j], V[i,j]])
    tablestrings.append(table.get_string())

print '\n'.join(tablestrings)

plt.clf()
plt.suptitle('Diagramm 5.1: Vergleich der Filterschaltungen')
barwidth = 0.35
ylim = [0, 1.5]
xticks = np.arange(len(f))
for i in range(len(labels)):
    ax = plt.subplot(2,3,i+1)
    plt.title(labels[i])
    plt.ylim(*ylim)
    p1_i = plt.bar(xticks, U[i], barwidth, color='b')
    p2_i = plt.bar(xticks+barwidth, V[i], barwidth, color='y')
    if i==0:
        p1 = p1_i
        p2 = p2_i
    ax.set_xticks(xticks+barwidth)
    ax.set_xticklabels(['$f_1$', '$f_2$', '$f_3$'])
    plt.axhline(1, color='black')
plt.figlegend((p1, p2), (r'$U \, [\frac{V}{V_{rms}}]$ absolut', r'$\frac{U}{U_0}$ relativ'), 'upper left', prop={'size':10})
papstats.savefig_a4('5.1.png')


# Dämpfung des f_1-Signals beim RC-Hochpass

f_G = unc.ufloat(3.08, 0.05)*const.kilo
f_1 = unc.ufloat(103,6)

U_theo = 1./unc.umath.sqrt(1+(f_G/f_1)**2)

print papstats.pformat(U_theo, format='c')




