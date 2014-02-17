# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c
import scipy.optimize as opt
import scipy.integrate as int

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

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

sl_bounds = [450, 1e5]
sl = (f > sl_bounds[0]) & (f < sl_bounds[1])

def fit_gf(f, V, O1, O2, n1, n2):
    return V/np.sqrt(1+1/(f/O1)**(2*n1))/np.sqrt(1+(f/O2)**(2*n2))

popt, pcov = opt.curve_fit(fit_gf, f[sl], gf[sl], p0=[1000,1000,50000,5,5])
pstats = papstats.PAPStats(gf[sl], fit_gf(f[sl], *popt))

varstr = ['V','\Omega_1','\Omega_2','n_1','n_2']
varunit = [None,'Hz','Hz',None,None]

plt.clf()
plt.scatter(f, gf, label='Messpunkte')
fspace = np.logspace(np.log10(sl_bounds[0]),np.log10(sl_bounds[1]), num=200)
plt.plot(fspace, fit_gf(fspace, *popt), label=r'Fit $g(f)=\frac{V}{\sqrt{1+(\frac{\Omega_1}{f})^{2*n_1}}\sqrt{1+(\frac{f}{\Omega_2})^{2*n_2}}}$ mit:'+'\n%s\n%s' % (''.join(['\n'+papstats.pformat(popt[i],pcov[i,i],label=varstr[i],unit=varunit[i]) for i in range(len(popt))]), pstats.legendstring()))
#for sl_bound in sl_bounds:
#    plt.vlines(sl_bound, 5, 100, label='Fit-Bereich Grenze bei %s' % papstats.pformat(sl_bound, unit='Hz'))
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e1,1e6)
plt.ylim(4,2e3)
plt.legend(loc='lower center')
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('2.png', dpi=144)

# Berechnung der Bandbreite

def gf_sq(f, *p):
    return fit_gf(f, *p)**2

B = int.quad(gf_sq, 0, np.inf, args=tuple(popt))[0]
print "B = %f" % B


#####
print "3 (Rauschspannung)"
#####

data = np.loadtxt('2.1.txt', skiprows=1)

N = data[:,0]
R = data[:,1]
dRrel = 0.005 # max. 0.5% Fehler auf die Widerstände
dR = R * dRrel
Uaus = data[:,2]
dUrel_stat = 0.003 # 0.3% statistischer Fehler
dUaus = np.sqrt((Uaus*dUrel_stat)**2+(data[:,3]/np.sqrt(N))**2)  # statistischer Fehler Fehler des Mittelwerts TODO: quadratisch addieren?

Usq = Uaus[1:]**2-Uaus[0]**2
dUsq = np.sqrt((2*Uaus[1:]*dUaus[1:])**2+(2*Uaus[0]*dUaus[0])**2)

def fit_Usq(x, c):
    return c*x

popt, pcov = opt.curve_fit(fit_Usq, R[1:], Usq, sigma=np.sqrt(dR[1:]**2+dUsq**2))
c = popt[0]
dc = pcov[0,0]
pstats = papstats.PAPStats(Usq, fit_Usq(R[1:], c), sigma=np.sqrt(dR[1:]**2+dUsq**2))

plt.clf()
plt.errorbar(R[1:], Usq, xerr=dR[1:], yerr=dUsq, ls='none', label='Messpunkte')
plt.plot(R[1:], fit_Usq(R[1:], c), label='Fit $(U_aus^2-U_V^2)=c*R$ mit:\n%s\n%s' % (papstats.pformat(c, dc, unit='\frac{V}{\Omega}'), pstats.legendstring()))
fig = plt.gcf()
fig.set_size_inches(11.69,8.27)
plt.savefig('3.png', dpi=144)

