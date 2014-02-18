# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import scipy.optimize as opt
import scipy.integrate as int

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

# Raumtemperatur
T = 23+const.zero_Celsius
dT = 0.1

# Dämpfung
D = 0.001
dDrel = 0.2
dD = D * dDrel


#####
print "2 (Frequenzgang)"
#####

data = np.loadtxt('3.2.txt', skiprows=1)

Uein = 0.2
#dUein

f = data[:,0]
#df
Uaus = data[:,1]
#dUaus

gf = Uaus/Uein/D
#dgf

# Fit Frequenzgang

sl_bounds = [5e2, 1e5]
#sl_bounds = [1,1e7]

sl = (f > sl_bounds[0]) & (f < sl_bounds[1])

def fit_gf(f, V, O1, O2, n1, n2):
    return V/np.sqrt(1+1/(f/O1)**(2*n1))/np.sqrt(1+(f/O2)**(2*n2))

popt, pcov = opt.curve_fit(fit_gf, f[sl], gf[sl], p0=[1000,1000,50000,5,5])
pstats = papstats.PAPStats(gf[sl], fit_gf(f[sl], *popt), ddof=5)

varstr = ['V','\Omega_1','\Omega_2','n_1','n_2']
varunit = [None,'Hz','Hz',None,None]

# Berechnung der Bandbreite

def gf_sq(f, *p):
    return fit_gf(f, *p)**2

B = int.quad(gf_sq, 0, np.inf, args=tuple(popt))[0]
print B
#B = 48100807970.4

# plot
plt.clf()
plt.scatter(f, gf, label='Messpunkte', s=10, c='black', marker='s')
fspace = np.logspace(np.log10(sl_bounds[0]),np.log10(sl_bounds[1]), num=200)
plt.plot(fspace, fit_gf(fspace, *popt), label=r'Fit $g(f)=\frac{V}{\sqrt{1+(\frac{\Omega_1}{f})^{2*n_1}}\sqrt{1+(\frac{f}{\Omega_2})^{2*n_2}}}$ mit:'+'\n%s\n%s\n$B=\int_0^\infty \! g(f)^2 \mathrm{d}f=$%s' % (''.join(['\n'+papstats.pformat(popt[i],pcov[i,i],label=varstr[i],unit=varunit[i]) for i in range(len(popt))]), pstats.legendstring(), papstats.pformat(B, signi=3, unit='Hz')))
ybounds = (5,1.5e3)
for i in range(len(sl_bounds)):
    plt.annotate('Fit-Bereich Grenze\nbei %s' % papstats.pformat(sl_bounds[i],  unit='Hz'), (sl_bounds[i], ybounds[1]), ha=['right','left'][i], va='top', textcoords='offset points',  xytext=([-1,1][i]*10,0))
    plt.vlines(sl_bounds[i], *ybounds)
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e1,1e6)
plt.ylim(4,2e3)
plt.legend(loc='lower center')
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('2.png', dpi=144)


#####
print "3 (Boltzmannkonstante)"
#####

# Messdaten der Rauschspannungsmessung

data = np.loadtxt('2.1.txt', skiprows=1)

N = data[:,0]
R = data[:,1]*const.kilo
dRrel = 0.005 # max. 0.5% Fehler auf die Widerstände
dR = R * dRrel
Uaus = data[:,2]*const.milli
dUrel_stat = 0.003 # 0.3% statistischer Fehler
dUaus = np.sqrt((Uaus*dUrel_stat)**2+(data[:,3]*const.milli/np.sqrt(N))**2)  # statistischer Fehler Fehler des Mittelwerts TODO: quadratisch addieren?

Usq = Uaus[1:]**2-Uaus[0]**2
dUsq = np.sqrt((2*Uaus[1:]*dUaus[1:])**2+(2*Uaus[0]*dUaus[0])**2)

# Linearer Fit Usq=c*R

def fit_Usq(x, c):
    return c*x

popt, pcov = opt.curve_fit(fit_Usq, R[1:], Usq, sigma=np.sqrt(dR[1:]**2+dUsq**2))
c = popt[0]
dc = pcov[0,0]
pstats = papstats.PAPStats(Usq, fit_Usq(R[1:], c), sigma=np.sqrt(dR[1:]**2+dUsq**2), ddof=1)
print pstats

plt.clf()
plt.errorbar(R[1:], Usq, xerr=dR[1:], yerr=dUsq, ls='none', label='Messpunkte')
plt.plot(R[1:], fit_Usq(R[1:], c), label='Fit $(U_{aus}^2-U_V^2)=c*R$ mit:\n%s\n%s' % (papstats.pformat(c, dc, label='c', unit=r'\frac{V}{\Omega}'), pstats.legendstring()))
plt.legend()
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.png', dpi=144)

# Berechnung der Boltzmannkonstante

kB = c/(4*T*B)
dkBrel_stat = np.sqrt((dc/c)**2+(dT/T)**2)
dkB_stat = kB * dkBrel_stat
dkBrel_syst = 0.02
dkB_syst = kB * dkBrel_syst
print np.array([kB, dkB_stat, dkB_syst])*1e23