# -*- coding: utf-8 -*-

__author__ = 'christian'

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c
import scipy.optimize as opt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

def setAxesLabel():
    plt.xlabel('U_E in V')
    plt.ylabel('U_A in V')

ohm = u"\u03A9"

plt.ion()
plt.cla()

# Eingangswiderstand
R_E = 3 * c.kilo
dR_rel = 5e-2
dR_E = R_E * dR_rel

# Gegenkopplungswiderstand
R_G = np.array([48.7, 274, 680]) * c.kilo
dR_G = dR_rel * R_G

print "Unsere Widerstände:"
print "Eingangswiderstand:", R_E, '+/-', dR_E
print "Gegenkopplungswiderstände:", R_G, "+/-", dR_G

print "Messung 1 a)"

# Laden der Messdaten
data = np.loadtxt('1.a.txt')

U_E = data[:,0]
dU_E = 2 * c.milli

U_A = data[:, 1:3]
dU_A_rel = 3e-2
dU_A = U_A * dU_A_rel + 40*c.milli

# Plotten der Messwerte mit den Fehlern
plt.errorbar(x=U_E, y=U_A[:,0], xerr=dU_E, yerr=dU_A[:,0],ls='none', label='R_G=48.7kOhm')
plt.errorbar(x=U_E, y=U_A[:,1], xerr=dU_E, yerr=dU_A[:,1],ls='none', label='R_G=274kOhm')
plt.title('Messdaten mit Fehlern')
setAxesLabel()
plt.legend()

plt.savefig('1.a.messdaten.png')

plt.cla()


# Fitting der Daten:
# V_0 = - U_A / U_E
# => U_A = -V_0 * U_E
V_0, dV_0 = (np.zeros(2), np.zeros(2))

pktSelect = [slice(None, None), slice(2, -3)]

def fU_A(U_E, V_0):
    return -V_0 * U_E


# Fitting für 48.7kOhm
popt, pcov = opt.curve_fit(fU_A, U_E, U_A[pktSelect[0],0], sigma=dU_A[pktSelect[0],0])
V_0[0]  = popt[0]
dV_0[0] = pcov[0,0]

# Fitting für 247kOhm
popt, pcov = opt.curve_fit(fU_A, U_E[pktSelect[1]], U_A[pktSelect[1],1], sigma=dU_A[pktSelect[1],1])
V_0[1]  = popt[0]
dV_0[1] = pcov[0,0]

print "Fit-Werte:"
print "V_0:", V_0, "+/-", dV_0


# Berechnete Verstärkung
V_b = R_G[0:2] / R_E
dV_b_rel = np.sqrt(2) * dR_rel
dV_b = dV_b_rel * V_b

print "Berechnete Betriebsverstärkung:"
print "V_b:", V_b, "+/-", dV_b

# ChiSquare

pstats = [papstats.PAPStats(U_A[pktSelect[i],i], fU_A(U_E[pktSelect[i]], V_0[i]), dU_A[pktSelect[i], i], 1) for i in range(2)]

# Plot 1 48.7kOhm
plt.errorbar(x=U_E, y=U_A[:,0], xerr=dU_E, yerr=dU_A[:,0],ls='none')
plt.plot(U_E, fU_A(U_E, V_0[0]))
plt.title(ur'Ausgangsspannung als Funktion der Eingangsspannung für R_G = 48.7k$\Omega$')
setAxesLabel()

plt.savefig('1.a.1.png')

plt.cla()


# Plot 2 274kOhm
plt.errorbar(x=U_E, y=U_A[:,1], xerr=dU_E, yerr=dU_A[:,1],ls='none')
plt.plot(U_E[pktSelect[1]], fU_A(U_E[pktSelect[1]], V_0[1]))
plt.title(ur'Ausgangsspannung als Funktion der Eingangsspannung für R_G = 274k$\Omega$')
setAxesLabel()

plt.savefig('1.a.2.png')