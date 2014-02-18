# -*- coding: utf-8 -*-

import numpy as np
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as const
import scipy.optimize as opt
import scipy.integrate as int

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

# Raumtemperatur
T = unc.ufloat(23+const.zero_Celsius, 0.1, 'T')

# Dämpfung
D = 0.001

#####
print "2 (Frequenzgang)"
#####

data = np.loadtxt('3.2.txt', skiprows=1)

Uein = 0.2
f = data[:,0]
Uaus = data[:,1]
Uaus -= 0.0014073 # Verstärkerrauschen

gf = Uaus/Uein/D

# Fit Frequenzgang

def fit_gf(f, V, O1, O2, n1, n2):
    return V/np.sqrt(1+1/(f/O1)**(2*n1))/np.sqrt(1+(f/O2)**(2*n2))

# Optimalen Fit-Bereich finden
'''
sl_bounds = [1,1e7]
opt_chisq = None
opt_bounds = [1,1e7]

for i in range(2):
    if i is 0:
        sl_space = np.logspace(2, 3, num=100)
    else:
        sl_space = np.logspace(4.7, 5.5, num=100)
    for sl_bound in sl_space:
        sl_bounds[i] = sl_bound
        
        sl = (f > sl_bounds[0]) & (f < sl_bounds[1])

        popt, pcov = opt.curve_fit(fit_gf, f[sl], gf[sl], p0=[1000,1000,50000,5,5])
        pstats = papstats.PAPStats(gf[sl], fit_gf(f[sl], *popt), ddof=5)

        if opt_chisq is None or pstats.chisq[1] < opt_chisq:
            opt_chisq = pstats.chisq[1]
            opt_bounds[i] = sl_bound
            print opt_chisq, opt_bounds

    opt_chisq = None

sl_bounds = opt_bounds
'''

sl_bounds = [4.5e2, 1e5] # passender Fit-Bereich
sl = (f > sl_bounds[0]) & (f < sl_bounds[1])

# Fit mit Startwerten
popt, pcov = opt.curve_fit(fit_gf, f[sl], unp.nominal_values(gf[sl]), p0=[1000,1000,50000,5,5])
pstats = papstats.PAPStats(unp.nominal_values(gf[sl]), fit_gf(f[sl], *popt), ddof=5)

popt = [unc.ufloat(popt[i], pcov[i,i]) for i in range(len(popt))]
plabel = ['V','\Omega_1','\Omega_2','n_1','n_2']
punit = [None,'Hz','Hz',None,None]

# Berechnung der Bandbreite

def gf_sq(f, *p):
    return fit_gf(f, *p)**2

# Integral der Quadratfunktion von g(f) von 0 bis Unendlich
B = int.quad(gf_sq, 0, np.inf, args=tuple(unp.nominal_values(popt)))[0]
B = unc.ufloat(B, B*0.02, 'B') # 2% systematischer Fehler auf B

# plot
plt.clf()
#plt.title(u'Diagramm 3.1: Frequenzgang der Messelektronik mit Verstärker und Bandfilter')
plt.title(u'Diagramm 3.2: Frequenzgang der Messelektronik mit Verstärker und Bandfilter (Korrektur)')
plt.scatter(f, gf, label='Messpunkte', s=10, c='black', marker='s')
fspace = np.logspace(np.log10(sl_bounds[0]),np.log10(sl_bounds[1]), num=200)
plt.plot(fspace, fit_gf(fspace, *unp.nominal_values(popt)), label=r'Fit $g(f)=\frac{V}{\sqrt{1+(\frac{\Omega_1}{f})^{2*n_1}}\sqrt{1+(\frac{f}{\Omega_2})^{2*n_2}}}$ mit:'+'\n%s\n%s\n$B=\int_0^\infty \! g(f)^2 \mathrm{d}f=$%s\n\n%s' % (''.join(['\n'+papstats.pformat(popt[i], label=plabel[i], unit=punit[i]) for i in range(len(popt))]), pstats.legendstring(), papstats.pformat(B, unit='Hz'), r'(mit $2\%$ system. Fehler)'))
ybounds = (5,1.5e3)
for i in range(len(sl_bounds)):
    plt.annotate('Fit-Bereich Grenze\nbei %s' % papstats.pformat(sl_bounds[i],  unit='Hz'), (sl_bounds[i], ybounds[1]), ha=['right','left'][i], va='top', textcoords='offset points',  xytext=([-1,1][i]*10,0))
    plt.vlines(sl_bounds[i], *ybounds)
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e1,1e6)
plt.ylim(2.5,2e3)
plt.xlabel('Frequenz $f \, [Hz]$')
plt.ylabel('$g(f)$')
plt.legend(loc='lower center')
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.png', dpi=144)


#####
print "3 (Boltzmannkonstante)"
#####

# Messdaten der Rauschspannungsmessung

data = np.loadtxt('2.1.txt', skiprows=1)

N = data[:,0]
R = data[:,1]*const.kilo
dRrel = 0.005 # max. 0.5% Fehler auf die Widerstände
R = unp.uarray(R, R*dRrel)
Uaus = data[:,2]*const.milli
dUrel_mess = unc.ufloat(1, 0.003) # 0.3% Messfehler
dUrel_stat = unp.uarray(1, (data[:,3]*const.milli/np.sqrt(N)) / Uaus) # Fehler des Mittelwerts
Uaus = unp.uarray(Uaus, Uaus * unp.std_devs(dUrel_mess + dUrel_stat))  # Messfehler und Fehler des Mittelwerts quadratisch addiert

Usq = Uaus[1:]**2-Uaus[0]**2

# Linearer Fit Usq=c*R

def fit_Usq(x, c):
    return c*x

popt, pcov = opt.curve_fit(fit_Usq, unp.nominal_values(R[1:]), unp.nominal_values(Usq), sigma=unp.std_devs(R[1:]+Usq))
c = unc.ufloat(popt[0], pcov[0,0], 'c')
pstats = papstats.PAPStats(unp.nominal_values(Usq), fit_Usq(unp.nominal_values(R[1:]), c.n), sigma=unp.std_devs(R[1:]+Usq), ddof=1)

plt.clf()
plt.title('Diagramm 3.3: Nyquist Beziehung zwischen effektiver Rauschspannung und Widerstand')
plt.errorbar(unp.nominal_values(R[1:])/const.kilo, unp.nominal_values(Usq)/(const.milli**2), xerr=unp.std_devs(R[1:])/const.kilo, yerr=unp.std_devs(Usq)/(const.milli**2), ls='none', label='Messpunkte')
plt.plot(unp.nominal_values(R[1:])/const.kilo, fit_Usq(unp.nominal_values(R[1:]), c.n)/(const.milli**2), label='Fit $(U_{aus}^2-U_V^2)=c*R$ mit:\n%s\n%s' % (papstats.pformat(c, label='c', unit=r'\frac{V^2}{\Omega}'), pstats.legendstring()))
plt.legend(loc='lower right')
plt.xlabel('Widerstand $R \, [k\Omega]$')
plt.ylabel('$(U_{aus}^2-U_V^2) \, [mV^2]$')
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.3.png', dpi=144)

# Berechnung der Boltzmannkonstante

kB = c/(4*T*B)
print kB

papstats.print_rdiff(kB, unc.ufloat(1.3806488e-23, 1.3e-29))