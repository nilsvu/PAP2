import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c
import scipy.optimize as opt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstat

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

slice = [slice(None, None), slice(2, -2)]

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
chisq = np.zeros((2,2))
chisq[0] = np.array(papstat.chisquared(Ua[slice[0],0], lin(Ue[slice[0]], -V0[0]), std=dUa[slice[0],0], ddof=1))
chisq[1] = np.array(papstat.chisquared(Ua[slice[1],1], lin(Ue[slice[1]], -V0[1]), std=dUa[slice[1],1], ddof=1))

print chisq

# plot
plt.clf()
plt.errorbar(Ue, Ua[:,0], xerr=dUe, yerr=dUa[:,0], ls='none')
plt.plot(Ue, -V0[0] * Ue)
plt.errorbar(Ue, Ua[:,1], xerr=dUe, yerr=dUa[:,1], ls='none')
plt.plot(Ue, -V0[1] * Ue)
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('1.a.png', dpi=144)