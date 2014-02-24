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

popt, pstats = papstats.curve_fit(fit_Uind, w, Uind)

# Berechnung der Magnetfeldstärke
B = popt[0]/nF/AF
papstats.print_rdiff(B, B_helmh(I))

Uind_true = nF*B_helmh(I)*AF*w

plt.clf()
plt.title(u'Diagramm 3.1: Induktionsspannung in Abhängigkeit von der Rotationsfrequenz der Flachspule')
plt.xlabel('Kreisfrequenz $\omega \, [Hz]$')
plt.ylabel('Max. Induktionsspannung $U_{ind} \, [V]$')
papstats.plot_data(w, Uind, label='Messpunkte')
papstats.plot_fit(fit_Uind, popt, pstats, np.linspace(w[0].n, w[-1].n), eq='U_{ind}=c*\omega', punits=[r'T \times m^2'])
papstats.plot_data(w, Uind_true, label='Erwartungswerte '+r'$U_{ind}=\frac{8*\mu_0*n_H*I}{\sqrt{125}*R}*n_F*A_F*\omega$')
plt.legend(loc='lower right')
papstats.savefig_a4('3.1.png')

# Vergleich mit Erwartungswerten

V = Uind/Uind_true

sigma_upper = unp.nominal_values(V+unp.std_devs(V)*3)
sigma_lower = unp.nominal_values(V-unp.std_devs(V)*3)

plt.clf()
plt.title(u'Diagramm 3.1.b: Verhältnis aus gemessener und erwarteter Induktionsspannung mit $3 \sigma$-Bereich')
plt.xlabel('Kreisfrequenz $\omega \, [Hz]$')
plt.ylabel(u'Verhältnis '+r'$\frac{U_{ind}^{mess}}{U_{ind}^{erw}}$')
papstats.plot_data(w, V, label='Messpunkte')
plt.fill_between(unp.nominal_values(w), sigma_upper, sigma_lower, color='red', alpha=0.5)
plt.axhline(y=1., c='black')
plt.ylim(0.5, 1.5)
papstats.savefig_a4('3.1.b.png')

# b) Abhängigkeit vom Spulenstrom

w = unc.ufloat(10.1, 0.1)*2*const.pi

data = np.loadtxt('2.b.txt', skiprows=1)

I = unp.uarray(data[:,0], data[:,1])
Uind = unp.uarray(data[:,2], data[:,3])/2

popt, pstats = papstats.curve_fit(fit_Uind, I, Uind)

Uind_true = nF*B_helmh(I)*AF*w

plt.clf()
plt.title(u'Diagramm 3.2: Induktionsspannung in Abhängigkeit vom Spulenstrom')
plt.xlabel('Spulenstrom $I \, [A]$')
plt.ylabel('Max. Induktionsspannung $U_{ind} \, [V]$')
papstats.plot_data(I, Uind, label='Messpunkte')
papstats.plot_fit(fit_Uind, popt, pstats, np.linspace(I[0].n, I[-1].n), eq='U_{ind}=c*I', punits=[r'\frac{T \times m^2}{A \times s}'])
papstats.plot_data(I, Uind_true, label='Erwartungswerte '+r'$U_{ind}=\frac{8*\mu_0*n_H*I}{\sqrt{125}*R}*n_F*A_F*\omega$')
plt.legend(loc='lower right')
papstats.savefig_a4('3.2.png')


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

popt, pstats = papstats.curve_fit(fit_cos, a, Uind)

plt.clf()
plt.subplot(121)
plt.xlabel(r'Winkel $\alpha \, [\pi]$')
plt.ylabel('Max. Induktionsspannung $U_{ind} \, [V]$')
papstats.plot_data(a/const.pi, Uind, label='Messpunkte')
papstats.plot_fit(fit_cos, popt, pstats, np.linspace(0,2*const.pi,num=100), xscale=1./const.pi, eq=r'U_{ind}=c*|cos(\alpha + \phi)|', plabels=['c',r'\phi'], punits=['V','rad'])
plt.xlim(0,2)
plt.ylim(0,1.6)
plt.legend(loc='upper center')
# Kreisplot
plt.subplot(122, polar=True)
papstats.plot_data(a, Uind, capsize=0)
papstats.plot_fit(fit_cos, popt, pstats, np.linspace(0,2*const.pi,num=100))
plt.xlabel(r'Winkel $\alpha$')
plt.ylabel('Max. Induktionsspannung $U_{ind} \, [V]$')
fig = plt.gcf()
fig.suptitle(u'Diagramm 3.3: Induktionsspannung in Abhängigkeit vom Winkel der Flachspule')
papstats.savefig_a4('3.3.png')

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
papstats.plot_data(W/const.kilo, V, label='Messpunkte')
plt.xlabel('Wechselstromfrequenz '+r'$\Omega \, [kHz]$')
plt.ylabel(u'Verhältnis '+r'$\frac{U_{ind}}{U_H}$')
plt.legend()
papstats.savefig_a4('3.4.png')

# c) Widerstand

R = Uhelmh/Ihelm

def fit_R(W, L):
    return L*W

popt, pstats = papstats.curve_fit(fit_R, W, R)

plt.clf()
plt.title(u'Diagramm 3.5: Widerstand der Helmholtzspulen als Funktion der Frequenz')
papstats.plot_data(W/const.kilo, R, label='Messpunkte')
papstats.plot_fit(fit_R, popt, pstats, np.linspace(W[0].n, W[-1].n), xscale=1./const.kilo, eq=r'Z=L*\Omega', punits=['H'])
plt.xlabel('Wechselstromfrequenz '+r'$\Omega \, [kHz]$')
plt.ylabel(u'Widerstand '+r'$Z=\frac{U_H}{I_H} \, [\Omega]$')
plt.legend(loc='lower right')
papstats.savefig_a4('3.5.png')


#####
print('\n# 3 Bestimmung des Erdmagnetfeldes durch Kompensation')
#####

Uind = unc.ufloat(170,15)*const.milli/2
w = unc.ufloat(14.5,0.1)*2*const.pi

Berd1 = Uind/(nF*AF*w)
papstats.print_rdiff(Berd1/const.micro, unc.ufloat(49,0))

# Kompensationsmessung
Uind = unc.ufloat(94.4,5)*const.milli/2
Ihelm = unc.ufloat(62.8,0.05)*const.milli
Bhor = Uind/(nF*AF*w)
Bver = B_helmh(Ihelm)
print Bhor, Bver
Berd2 = unc.umath.sqrt((Bver**2)+(Bhor**2))
papstats.print_rdiff(Berd2/const.micro, unc.ufloat(49,0))
papstats.print_rdiff(Berd1/const.micro, Berd2/const.micro)

a = unc.umath.atan(Bver/Bhor) # Inklinationswinkel
papstats.print_rdiff(a/(2*const.pi)*360, unc.ufloat(66,0))
