import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c
import scipy.optimize as opt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

# Eingangswiderstand
Re = 3*c.kilo
dRrel = 0.05
dRe = Re * dRrel

# Gegenkopplungswiderstaende
Rg = np.array([48.7, 274, 680])*c.kilo
dRg = dRrel * Rg


#####
print "# 1.a"
#####

# Spannungsmessdaten

data = np.loadtxt('1.a.txt')

Ue = data[:,0]
dUe = 2*c.milli
dUrel = 0.03
Ua = data[:,1:3]
dUa = dUrel*Ua + 2*c.milli

# Linearer Fit
# V0 = -Ua/Ue
# => Ua = -V0 * Ue
# => V0 ist negative Steigung der Fit-Geraden

V0 = np.zeros(2)
dV0 = np.zeros(2)

slice = [slice(None, None), slice(2, -3)]

def lin(x, m):
    return m * x

popt, pcov = opt.curve_fit(lin, Ue[slice[0]], Ua[slice[0],0], sigma=dUa[slice[0],0])
V0[0] = popt[0]
dV0[0] = pcov[0,0]
popt, pcov = opt.curve_fit(lin, Ue[slice[1]], Ua[slice[1],1], sigma=dUa[slice[1],1])
V0[1] = popt[0]
dV0[1] = pcov[0,0]

V0 *= -1

print "Linearer Fit Ua = -V0 * Ue:"
print "V0 = "+str(V0)+"+/-"+str(dV0)

# Berechnete Betriebsverstaerkung
Vb = Rg[0:2]/Re
dVbrel = np.sqrt(2)*dRrel
dVb = dVbrel * Vb

print "Berechnete Betriebsverstaerkung:"
print "Vb = "+str(Vb)+"-/+"+str(dVb)

# Chi-Quadrate
pstats = [papstats.PAPStats(Ua[slice[i],i], lin(Ue[slice[i]], -V0[i]), std=dUa[slice[i],i], ddof=1) for i in range(2)]

print "PAPStats:", pstats

# plot
plt.clf()
plt.title('Ausgangsspannung als Funktion der Eingangsspannung fuer Gleichstrom')
plt.errorbar(Ue, Ua[:,0], xerr=dUe, yerr=dUa[:,0], ls='none', label=r'Rg=Rg1=48.7k$\Omega$')
plt.plot(Ue, -V0[0] * Ue, label='Fit zu Rg=Rg1')
plt.errorbar(Ue, Ua[:,1], xerr=dUe, yerr=dUa[:,1], ls='none', label=r'Rg=Rg2=274k$\Omega$')
plt.plot(Ue, -V0[1] * Ue, label='Fit zu Rg=Rg2')
plt.xlabel('Ue (Eingangsspannung) [V]')
plt.ylabel('Ua (Ausgangsspannung) [V]')
plt.legend()
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('1.a.png', dpi=144)