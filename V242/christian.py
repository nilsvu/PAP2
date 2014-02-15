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

plt.ion()
plt.cla()

fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
fig.set_dpi(80)

# Eingangswiderstand
R_E = 3 * c.kilo
dR_rel = 5e-2
dR_E = R_E * dR_rel

# Gegenkopplungswiderstand
R_G = np.array([48.7, 274, 680]) * c.kilo
dR_G = dR_rel * R_G

print "Unsere Widerstände:"
print "Eingangswiderstand:", R_E, '+/-', dR_E
print "Gegenkopplungswiderstände:"
for i in range(R_G.size):
    print R_G[i], "+/-", dR_G[i]


print ''

print "Messung 1 a)"

# Laden der Messdaten
data = np.loadtxt('1.a.txt')

U_E = data[:,0]
dU_E = 2 * c.milli

U_A = data[:, 1:3]
dU_A_rel = 3e-2
dU_A = U_A * dU_A_rel + 40*c.milli

# Plotten der Messwerte mit den Fehlern
plt.errorbar(x=U_E, y=U_A[:,0], xerr=dU_E, yerr=dU_A[:,0],ls='none', label=r'R_G=48.7k$\Omega$')
plt.errorbar(x=U_E, y=U_A[:,1], xerr=dU_E, yerr=dU_A[:,1],ls='none', label=r'R_G=274k$\Omega$')
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

for pstat in pstats:
    print pstat

plt.cla()





################################################################
################################################################
##############          1 b)                ####################
################################################################
################################################################

print ''
print 'Messung 1b'

data = np.loadtxt('1.b.txt')

U_E = data[:,0] / 10.
dU_E_rel = 3e-2
dU_E = dU_E_rel * U_E

U_A = data[:, 1:3]
dU_A_rel = 3e-2
dU_A = U_A * dU_A_rel + 40*c.milli

plt.errorbar(U_E, U_A[:,0], xerr=dU_E, yerr=dU_A[:,0], ls='none', label=r'R_G=274k$\Omega$')
plt.errorbar(U_E, U_A[:,1], xerr=dU_E, yerr=dU_A[:,1], ls='none', label=r'R_G=680k$\Omega$')
setAxesLabel()
plt.legend(loc=2)

plt.savefig('1.b.messdaten.png')
plt.cla()

# Fitting der Daten:
# V_0 = - U_A / U_E
# => U_A = -V_0 * U_E

def fU_A(U_E, V_0):
    return V_0 * U_E

V_0, dV_0 = (np.zeros(2), np.zeros(2))

pktSelect = [slice(None, None), slice(2, None)]

# Fit für 274kOhm
popt, pcov = opt.curve_fit(fU_A, U_E[pktSelect[0]], U_A[pktSelect[0],0], sigma=dU_A[pktSelect[0],0])
V_0[0]  = popt[0]
dV_0[0] = pcov[0,0]

# Fit für 680kOhm
popt, pcov = opt.curve_fit(fU_A, U_E[pktSelect[1]], U_A[pktSelect[1],1], sigma=dU_A[pktSelect[1],1])
V_0[1]  = popt[0]
dV_0[1] = pcov[0,0]

print "Fit-Werte:"
print "V_0:", V_0, "+/-", dV_0

# Berechnete Verstärkung
V_b = R_G[1:3] / R_E
dV_b_rel = np.sqrt(2) * dR_rel
dV_b = dV_b_rel * V_b

print "Berechnete Betriebsverstärkung:"
print "V_b:", V_b, "+/-", dV_b

pstats = [papstats.PAPStats(U_A[pktSelect[i],i], fU_A(U_E[pktSelect[i]], V_0[i]), dU_A[pktSelect[i], i], 1) for i in range(2)]

# Plot 1 48.7kOhm
plt.errorbar(x=U_E, y=U_A[:,0], xerr=dU_E, yerr=dU_A[:,0],ls='none')
plt.plot(U_E[pktSelect[0]], fU_A(U_E[pktSelect[0]], V_0[0]))
plt.title(ur'Ausgangsspannung als Funktion der Eingangsspannung für R_G = 274k$\Omega$')
setAxesLabel()

plt.savefig('1.b.1.png')

plt.cla()


# Plot 2 274kOhm
plt.errorbar(x=U_E, y=U_A[:,1], xerr=dU_E, yerr=dU_A[:,1],ls='none')
plt.plot(U_E[pktSelect[1]], fU_A(U_E[pktSelect[1]], V_0[1]))
plt.title(ur'Ausgangsspannung als Funktion der Eingangsspannung für R_G = 680k$\Omega$')
setAxesLabel()

plt.savefig('1.b.2.png')
plt.cla()

for pstat in pstats:
    print pstat


################################################################
################################################################
##############          2                   ####################
################################################################
################################################################

data = np.loadtxt('2.a.txt')

f = data[:,0]
df = data[:, 1]

U_E = np.array([1, 0.3, 0.3]) / 10.
dU_E = U_E * 3e-2

U_A = data[:,2:5]
dU_A = U_A * 3e-2

V = U_A / U_E
dV = V * np.sqrt(2) * 3e-2

plt.xscale('log')
plt.yscale('log')
for i in range(V.shape[1]):
    plt.errorbar(f, V[:,i], xerr=df, yerr=dV[:,i], label='R_G='+str(R_G[i]/c.kilo)+r'k$\Omega$')
    pass

plt.xlabel('Frequenz in Hz')
plt.ylabel('V')



data = np.loadtxt('2.b.txt')

f   = data[:,0]
df = data[:,1]

U_E = 1. / 10.
dU_E = U_E * 3e-2

U_A = data[:,2]
dU_A = U_A * 3e-2

V = U_A / U_E
dV = V * np.sqrt(2) * 3e-2

plt.errorbar(f, V, xerr=df, yerr=dV, label=ur"R_G=48.7k$\Omega$ und C=560pF (rückgekoppelt)")

data = np.loadtxt('2.c.txt')

f   = data[:,0]
df = data[:,1]

U_E = 1. / 10.
dU_E = U_E * 3e-2

U_A = data[:,2]
dU_A = U_A * 3e-2

V = U_A / U_E
dV = V * np.sqrt(2) * 3e-2

plt.errorbar(f, V, xerr=df, yerr=dV, label=ur"R_G=48.7k$\Omega$ und C=47nF (vorgeschaltet)")


plt.legend(loc=1)
plt.savefig('2.png')