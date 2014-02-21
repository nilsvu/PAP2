# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 245 - Induktion
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as const
import scipy.optimize as opt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats


# Daten der Helmholtzspule

RH = 295*const.milli/2
aH = 147*const.milli
nH = 124

def B_helmh(I):
    return 8./np.sqrt(125)*const.mu_0*nH*I/RH

# Daten der Flachspule

AF = 41.7*const.centi**2
nF = 4e3

#####
print('# 1 (Induktionsgesetz)')
#####

# Fit Uind = n*B*A*w

def fit_Uind(x, c):
    return c*x

# a) Abhängigkeit von der Frequenz

I = unc.ufloat(4, 0.05)

data = np.loadtxt('2.a.txt', skiprows=1)

w = unp.uarray(data[:,0], data[:,1])*2*const.pi
Uind = unp.uarray(data[:,2], data[:,3])/2

popt, pcov = opt.curve_fit(fit_Uind, unp.nominal_values(w), unp.nominal_values(Uind), sigma=unp.std_devs(w+Uind))
popt = unp.uarray(popt, np.sqrt(np.diagonal(pcov)))
pstats = papstats.PAPStats(unp.nominal_values(w), fit_Uind(unp.nominal_values(w), *unp.nominal_values(popt)), sigma=unp.std_devs(w+Uind), ddof=1)

# Berechnung der Magnetfeldstärke
# mit c = n*B*A
B = popt[0]/nF/AF
print 'B =', B/const.milli, 'mT'
print 'Erwartungswert:\nB =', B_helmh(I)/const.milli, 'mT'
papstats.print_rdiff(B, B_helmh(I))

Uind_true = nF*B_helmh(I)*AF*w

plt.clf()
plt.title(u'Diagramm 3.1: Induktionsspannung in Abhängigkeit von der Rotationsfrequenz der Flachspule')
plt.xlabel('Kreisfrequenz $\omega \, [Hz]$')
plt.ylabel('Max. Induktionsspannung $U_{ind} \, [V]$')
plt.errorbar(unp.nominal_values(w), unp.nominal_values(Uind), xerr=unp.std_devs(w), yerr=unp.std_devs(Uind), ls='none', label='Messpunkte')
plt.plot(unp.nominal_values(w), fit_Uind(unp.nominal_values(w), *unp.nominal_values(popt)), label='Fit $U_{ind}=c*\omega$ mit:\n%s\n%s' % (papstats.pformat(popt[0], label='c'), pstats.legendstring()))
plt.errorbar(unp.nominal_values(w), unp.nominal_values(Uind_true), xerr=unp.std_devs(w), yerr=unp.std_devs(Uind_true), label='Erwartungswerte', ls='none')
plt.legend()
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.1.png', dpi=144)

# Vergleich mit Erwartungswerten

V = Uind/Uind_true

sigma_upper = unp.nominal_values(V+unp.std_devs(V)*3)
sigma_lower = unp.nominal_values(V-unp.std_devs(V)*3)

plt.clf()
plt.title(u'Diagramm 3.1.b: Verhältnis aus gemessener und erwarteter Induktionsspannung mit $3 \sigma$-Bereich')
plt.xlabel('Kreisfrequenz $\omega \, [Hz]$')
plt.ylabel(u'Verhältnis '+r'$\frac{U_{ind}^{mess}}{U_{ind}^{erw}}$')
plt.errorbar(unp.nominal_values(w), unp.nominal_values(V), xerr=unp.std_devs(w), yerr=unp.std_devs(V), label='Messpunkte')
plt.plot(unp.nominal_values(w), sigma_upper)
plt.plot(unp.nominal_values(w), sigma_lower)
plt.axhline(y=1., c='black')
plt.ylim(0.5, 1.5)
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.1.b.png', dpi=144)


# b) Abhängigkeit vom Spulenstrom

w = unc.ufloat(10.1, 0.1)*2*const.pi

data = np.loadtxt('2.b.txt', skiprows=1)

I = unp.uarray(data[:,0], data[:,1])
Uind = unp.uarray(data[:,2], data[:,3])/2

popt, pcov = opt.curve_fit(fit_Uind, unp.nominal_values(I), unp.nominal_values(Uind), sigma=unp.std_devs(I+Uind))
popt = unp.uarray(popt, np.diagonal(pcov))
pstats = papstats.PAPStats(unp.nominal_values(I), fit_Uind(unp.nominal_values(I), *unp.nominal_values(popt)), sigma=unp.std_devs(I+Uind), ddof=1)

Uind_true = nF*B_helmh(I)*AF*w

plt.clf()
plt.title(u'Diagramm 3.2: Induktionsspannung in Abhängigkeit vom Spulenstrom')
plt.xlabel('Spulenstrom $I \, [A]$')
plt.ylabel('Max. Induktionsspannung $U_{ind} \, [V]$')
plt.errorbar(unp.nominal_values(I), unp.nominal_values(Uind), xerr=unp.std_devs(I), yerr=unp.std_devs(Uind), ls='none', label='Messpunkte')
plt.plot(unp.nominal_values(I), fit_Uind(unp.nominal_values(I), *unp.nominal_values(popt)), label='Fit $U_{ind}=c*I$ mit:\n%s\n%s' % (papstats.pformat(popt[0], label='c'), pstats.legendstring()))
plt.errorbar(unp.nominal_values(I), unp.nominal_values(Uind_true), xerr=unp.std_devs(I), yerr=unp.std_devs(Uind_true), label='Erwartungswerte', ls='none')
plt.legend()
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.2.png', dpi=144)


#####
print('\n# 2 Induktionsspannung bei periodischem Feldstrom')
#####

# a) Winkelabhängigkeit

W = unc.ufloat(104,1)*2*const.pi # Kreisfrequenz der Wechselspannung

data = np.loadtxt('3.a.txt', skiprows=1)

a = unp.uarray(data[:,0], 2)/360.*2*const.pi
Uind = unp.uarray(data[:,1], data[:,2])/2

def fit_cos(a, c, d):
    return c*np.abs(np.cos(a+d))

popt, pcov = opt.curve_fit(fit_cos, unp.nominal_values(a), unp.nominal_values(Uind), sigma=unp.std_devs(a+Uind))
popt = unp.uarray(popt, np.sqrt(np.diagonal(pcov)))
pstats = papstats.PAPStats(unp.nominal_values(Uind), fit_cos(unp.nominal_values(a), *unp.nominal_values(popt)), sigma=unp.std_devs(a+Uind), ddof=2)

plt.clf()
plt.title(u'Diagramm 3.3: Induktionsspannung in Abhängigkeit vom Winkel der Flachspule')
plt.xlabel(r'Winkel $\alpha \, [\pi]$')
plt.ylabel('Max. Induktionsspannung $U_{ind} \, [V]$')
plt.errorbar(unp.nominal_values(a)/const.pi, unp.nominal_values(Uind), unp.std_devs(a)/const.pi, unp.std_devs(Uind), label='Messpunkte', ls='none')
aspace = np.linspace(0,2*const.pi,num=100)
plt.plot(aspace/const.pi, fit_cos(aspace, *unp.nominal_values(popt)), label=r'Fit $U_{ind}=c*|cos(\alpha + \phi)|$ mit:'+'\n%s\n%s' % (papstats.pformat(popt[0], label='c'), pstats.legendstring()))
plt.axvline(x=0, c='black')
plt.axvline(x=2, c='black')
plt.xlim(-0.25, 2.25)
plt.legend()
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.3.png', dpi=144)

# b) induzierte und angelegte Spannung

a = unc.ufloat(0,2)

data = np.loadtxt('3.b.txt', skiprows=1)

W = unp.uarray(data[:,0], data[:,1])*2*const.pi
Uhelmh = unp.uarray(data[:,2], data[:,3])/2
Ihelm = unp.uarray(data[:,4], data[:,5])*const.milli
Uind = unp.uarray(data[:,6], 0.05)/2

V = Uind/Uhelmh

plt.clf()
plt.title(u'Diagramm 3.4: Verhältnis von induzierter und angelegter Spannung als Funktion der Frequenz')
plt.errorbar(unp.nominal_values(W), unp.nominal_values(V), xerr=unp.std_devs(W), yerr=unp.std_devs(V), label='Messpunkte')
plt.legend()
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.4.png', dpi=144)

# c) Widerstand

R = Uhelmh/Ihelm

def fit_R(W, L):
    return L*W

popt, pcov = opt.curve_fit(fit_R, unp.nominal_values(W), unp.nominal_values(R), sigma=unp.std_devs(W+R))
popt = unp.uarray(popt, np.sqrt(np.diagonal(pcov)))
pstats = papstats.PAPStats(unp.nominal_values(R), fit_R(unp.nominal_values(W), *unp.nominal_values(popt)), sigma=unp.std_devs(R+W), ddof=1)

plt.clf()
plt.title(u'Diagramm 3.5: Widerstand der Helmholtzspulen als Funktion der Frequenz')
plt.errorbar(unp.nominal_values(W), unp.nominal_values(R), xerr=unp.std_devs(W), yerr=unp.std_devs(R), label='Messpunkte', ls='none')
plt.plot(unp.nominal_values(W), unp.nominal_values(fit_R(W, *popt)), label='Fit\n%s' % (papstats.pformat(popt[0]/const.milli, label='L', unit='mH')))
plt.legend()
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.5.png', dpi=144)


#####
print('\n# 3 Bestimmung des Erdmagnetfeldes durch Kompensation')
#####

Uind = unc.ufloat(320,5)*const.milli/2
w = unc.ufloat(14.5,0.1)*2*const.pi

Berd = Uind/(nF*AF*w)
print 'B =', Berd/const.micro, 'muT'
papstats.print_rdiff(Berd/const.micro, unc.ufloat(49,0))

# Kompensationsmessung
Uind = unc.ufloat(94.4,2)*const.milli/2
Ihelm = unc.ufloat(3.48,0.01)*const.milli
Bhor = Uind/(nF*AF*w)
Bver = B_helmh(Ihelm)
Berd = unc.umath.sqrt((Bver**2)+(Bhor**2))
print 'B =', Berd/const.micro, 'muT'
papstats.print_rdiff(Berd/const.micro, unc.ufloat(49,0))
a = unc.umath.atan(Bver/Bhor)
print a/(2*const.pi)*360

