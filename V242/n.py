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
print "# 1.a (Gleichspannung)"
#####

# Spannungsmessdaten

data = np.loadtxt('1.a.txt')

Ue = data[:,0]
dUe = 2*c.milli
dUrel = 0.03
Ua = data[:,1:3]
dUa = dUrel*Ua + 40*c.milli

# Linearer Fit
# V0 = -Ua/Ue
# => Ua = -V0 * Ue
# => V0 ist negative Steigung der Fit-Geraden

V0 = np.zeros(2)
dV0 = np.zeros(2)

sl = [slice(None, None), slice(2, -3)]

def lin(x, m):
    return m * x

popt, pcov = opt.curve_fit(lin, Ue[sl[0]], Ua[sl[0],0], sigma=dUa[sl[0],0])
V0[0] = popt[0]
dV0[0] = pcov[0,0]
popt, pcov = opt.curve_fit(lin, Ue[sl[1]], Ua[sl[1],1], sigma=dUa[sl[1],1])
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
pstats = [papstats.PAPStats(Ua[sl[i],i], lin(Ue[sl[i]], -V0[i]), std=dUa[sl[i],i], ddof=1) for i in range(2)]

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

###########
print "# 1.b (Wechselspannung)"
###########

data = np.loadtxt('1.b.txt')

dUrel = 0.03
Ue = data[:,0]/10. # Spannung wird um Faktor 10 heruntergesetzt
dUe = Ue * dUrel
Ua = data[:,1:3]
dUa = dUrel * Ua + 40*c.milli

# Linearer Fit
# V0 = -Ua/Ue
# => Ua = -V0 * Ue
# => V0 ist negative Steigung der Fit-Geraden

V0 = np.zeros(2)
dV0 = np.zeros(2)

sl = [slice(None, None), slice(2, None)]

popt, pcov = opt.curve_fit(lin, Ue[sl[0]], Ua[sl[0],0], sigma=dUa[sl[0],0])
V0[0] = popt[0]
dV0[0] = pcov[0,0]
popt, pcov = opt.curve_fit(lin, Ue[sl[1]], Ua[sl[1],1], sigma=dUa[sl[1],1])
V0[1] = popt[0]
dV0[1] = pcov[0,0]

V0 *= -1

print "Linearer Fit Ua = -V0 * Ue:"
print "V0 = "+str(V0)+"+/-"+str(dV0)

# Berechnete Betriebsverstaerkung
Vb = Rg[1:3]/Re
dVbrel = np.sqrt(2)*dRrel
dVb = dVbrel * Vb

print "Berechnete Betriebsverstaerkung:"
print "Vb = "+str(Vb)+"-/+"+str(dVb)

# Chi-Quadrate
pstats = [papstats.PAPStats(Ua[sl[i],i], lin(Ue[sl[i]], -V0[i]), std=dUa[sl[i],i], ddof=1) for i in range(2)]

print "PAPStats:", pstats

# plot
plt.clf()
plt.title('Ausgangsspannung als Funktion der Eingangsspannung fuer Wechselstrom')
plt.errorbar(Ue, Ua[:,0], xerr=dUe, yerr=dUa[:,0], ls='none', label=r'Rg=Rg2=274k$\Omega$')
plt.plot(Ue, -V0[0] * Ue, label='Fit zu Rg=Rg2')
plt.errorbar(Ue, Ua[:,1], xerr=dUe, yerr=dUa[:,1], ls='none', label=r'Rg=Rg3=680k$\Omega$')
plt.plot(Ue, -V0[1] * Ue, label='Fit zu Rg=Rg3')
plt.xlabel('Ue (Eingangsspannung) [V]')
plt.ylabel('Ua (Ausgangsspannung) [V]')
plt.legend()
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('1.b.png', dpi=144)


###########
print "# 2 (Frequenzgang bei Gegenkopplung)"
###########

data_a = np.loadtxt('2.a.txt')
data_b = np.loadtxt('2.b.txt')
data_c = np.loadtxt('2.c.txt')

dUrel = 0.03
Ue_a = np.array([1., 0.3, 0.3])/10. # Spannung wird um Faktor 10 heruntergesetzt
dUe_a = Ue_a * dUrel
Ua_a = data_a[:,2:5]
dUa_a = dUrel * Ua_a
Ue_b = 1./10.
dUe_b = Ue_b * dUrel
Ua_b = data_b[:,2]
dUa_b = Ua_b * dUrel
Ue_c = Ue_b
dUe_c = Ue_c * dUrel
Ua_c = data_c[:,2]
dUa_c = Ua_c * dUrel

f_a = data_a[:,0]
df_a = data_a[:,1]
f_b = data_b[:,0]
df_b = data_b[:,1]
f_c = data_c[:,0]
df_c = data_c[:,1]

# Verstaerkung V = Ua/Ue
dVrel = np.sqrt(2)*dUrel
V_a = Ua_a / Ue_a
dV_a = V_a * dVrel
V_b = Ua_b / Ue_b
dV_b = V_b * dVrel
V_c = Ua_c / Ue_c
dV_c = V_c * dVrel

# plot
plt.clf()
plt.title('Verstaerkung in Abhaengigkeit von der Frequenz (Frequenzgang)')
plt.xscale('log')
plt.yscale('log')
for i in range(V_a.shape[1]):
    plt.errorbar(f_a, V_a[:,i], xerr=df_a, yerr=dV_a[:,i], label='Rg=Rg'+str(i+1)+'='+str(Rg[i]/c.kilo)+r'k$\Omega$')
plt.errorbar(f_b, V_b, xerr=df_b, yerr=dV_b, label=r'Rg=Rg1=48.7k$\Omega$, C=560pF gegengekoppelt')
plt.errorbar(f_c, V_c, xerr=df_c, yerr=dV_c, label=r'Rg=Rg1=48.7k$\Omega$, C=47nF vorgeschaltet')
plt.xlabel('Frequenz [Hz]')
plt.ylabel('Verstaerkung')
plt.legend()
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('2.a.png', dpi=144)
