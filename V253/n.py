# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 252 - Aktivierung von Indium und von Silber mit thermischen Neutronen
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import scipy as sp
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as const
import scipy.interpolate as ip
import datetime as dt
import prettytable as pt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

# Zählrohr
r_Z = 7*const.milli
l = 0.04


#####
print "# 1 Nulleffekt"
#####

nU = unc.ufloat(120,np.sqrt(120))/(5.*60.)
print nU


#####
print "\n# 2 Absorption von b-Strahlung in Aluminium"
#####

data = np.loadtxt('1.txt', skiprows=1)

x = unp.uarray(data[:,0]*const.milli, 0.05*const.milli)
N = data[:,2]
n = unp.uarray(N, np.sqrt(N))/data[:,1]
nUb = unc.ufloat(149,np.sqrt(149))/(5.*60.)

n = n-nUb

table = pt.PrettyTable()
table.add_column('x [mm]', x/const.milli, align='r')
table.add_column('(n-n_U^b) [1/s]', n, align='r')
with open('Resources/3.1.txt', 'w') as file:
    file.write(table.get_string())

# Abschätzung der Maximalreichweite aus dem Diagramm
x1 = 3.3*const.milli
x2 = 3.8*const.milli
xM = unc.ufloat(x1+np.abs(x1-x2)/2., np.abs(x1-x2)/2.)

# Berechnung der Maximalenergie
RbES = 0.13*const.gram/const.centi**2
rAl = 2.7*const.gram/const.centi**3
Rb = RbES + xM * rAl
print "Flächendichte:", papstats.pformat(Rb/(const.gram/const.centi**2), label='R^b', unit='g/cm^2', format='.2u')
print "Maximalenergie:"
EM = unc.ufloat(2.25,0.1) # Aus Diagramm in Versuchsanleitung abgelesen
papstats.print_rdiff(EM,unc.ufloat(2.274,0))

# Plot
plt.clf()
plt.title(u'Diagramm 3.1: '+r'$\beta$'+u'-Strahlung in Abhängigkeit der Absorberdicke mit Extrapolation')
plt.yscale('log', nonposy='clip')
plt.xlabel('Absorberdicke $x \, [mm]$')
plt.ylabel(u'korrigierte Zählrate '+r'$(n-n_U^{\beta}) \, [\frac{Ereignisse}{s}]$')
papstats.plot_data(x/const.milli, n, label='Messpunkte')
xspace = np.linspace(0,0.004,num=1000)
sl = (n > 0)
for k in range(1,6):
    extrapolater = ip.UnivariateSpline(unp.nominal_values(x[sl]), np.log(unp.nominal_values(n[sl])), w=1./np.log(unp.std_devs(x+n)[sl]), k=k)
    plt.plot(xspace/const.milli, np.exp(extrapolater(xspace)), label='$k='+str(k)+'$', alpha=0.3)
sl_upper = 10
sl = sl & (n <= sl_upper)
k = 5
extrapolater = ip.UnivariateSpline(unp.nominal_values(x[sl]), np.log(unp.nominal_values(n[sl])), k=k)
plt.plot(xspace/const.milli, np.exp(extrapolater(xspace)), label='$k='+str(k)+'$ mit $n \in (0,'+str(sl_upper)+']$')
plt.fill_betweenx(np.linspace(1e-3,1e2), x1/const.milli, x2/const.milli, color='g', alpha=0.2)
plt.xlim(xspace[0]/const.milli, xspace[-1]/const.milli)
plt.ylim(1e-3,1e2)
handles, labels = plt.gca().get_legend_handles_labels()
p = plt.Rectangle((0, 0), 1, 1, color='g', alpha=0.2)
handles.append(p)
labels.append(u'Abschätzung der Maximalreichweite:\n$'+papstats.pformat(xM/const.milli, label='x_M', unit='mm')+'$')
plt.legend(handles, labels, loc='lower left', title='Extrapolation mit Polynomen der Ordnung $k$:')
papstats.savefig_a4('3.1.png')


#####
print "\n# 3 Absorption von y-Strahlung in Blei"
#####

data = np.loadtxt('2.txt', skiprows=1)

x = data[:,0]*const.milli
N = data[:,1]
n = unp.uarray(N, np.sqrt(N))/60

n = n-nU

def fit_damp(x, mu, n_0):
    return n_0*np.exp(-mu*x)

popt, pstats = papstats.curve_fit(fit_damp, x, n)

rhoPb = 11.34*const.gram/const.centi**3
k = popt[0]/rhoPb
print "Materialunabhängiger Schwächungskoeffizient:", papstats.pformat(k/(const.centi**2/const.gram), label='k', unit='cm^2/g', format='.2u')
Ey = unc.ufloat(1.45, 0.15)
print "Energie der y-Quanten"
papstats.print_rdiff(Ey,unc.ufloat((1.173+1.333)/2,0))

plt.clf()
plt.title('Diagramm 3.2: '+r'$\gamma$'+u'-Strahlung in Abhängigkeit der Absorberdicke')
plt.yscale('log')
plt.xlabel('Absorberdicke $x \, [mm]$')
plt.ylabel(u'korrigierte Zählrate '+r'$(n-n_U) \, [\frac{Ereignisse}{s}]$')
ylim = [1.5,1e2]
xlim = np.array([0,55])*const.milli
plt.ylim(*ylim)
plt.xlim(*xlim/const.milli)
papstats.plot_data(x/const.milli, n, label='Messpunkte')
papstats.plot_fit(fit_damp, popt, pstats, np.linspace(*xlim), eq=r'(n-n_U)=n_0*e^{-\mu*x}', plabels=[r'\mu','n_0'], punits=[r'm^{-1}','s^{-1}'], xscale=1./const.milli)
plt.legend()
papstats.savefig_a4('3.2.png')


#####
print "\n# 4 Bestimmung der Aktivität"
#####

# Erwartungswert
A0 = unc.ufloat(3.7e6,0)
t = (dt.date(2014,2,27)-dt.date(2012,2,2)).days
Th = 5.27*365.242199
A_erw = A0*np.exp(-np.log(2)/Th*t)
print "Erwartungswert:", papstats.pformat(A_erw, label='A_erw', unit='s^-1', format='.2u')

# Messdaten
d, dd, N = np.loadtxt('3.txt', skiprows=1, unpack=True)
d = unp.uarray(d, dd)
d += unc.ufloat(60,0.5)-unc.ufloat(160,0.5)+unc.ufloat(119,0.5) # Messprozess des Präparatabstands
d = unp.uarray(unp.nominal_values(d), 5)
d *= const.milli
N = unp.uarray(N, np.sqrt(N))
n = N/2./60. # 2 Photonen pro Zerfallsereignis

# Berechnung der Aktivität und Vergleich mit Erwartungswerten

def vergleiche_aktivitaet(A, filename, e_soll=None, k=None, index=None):
    headlines = ['d [cm]', 'A'+('_k'+str(index) if index is not None else '')+' [1/s]', 'Abweichung [1/s]', 'rel. Abweichung', 'Sigma-Bereich']
    if k is not None:
        headlines.insert(2, 'Korrekturfaktor')
    if e_soll is not None:
        headlines.append('Epsilon-Soll')
    table = pt.PrettyTable(headlines)
    for i in range(len(A)):
        row = [d[i]/const.centi, A[i]]
        if k is not None:
            row.append(k[i])
        row.extend(papstats.rdiff(A[i], A_erw))
        if e_soll is not None:
            row.append(e_soll[i])
        table.add_row(row)
    with open('Resources/'+filename+'.txt', 'w') as file:
        file.write(table.get_string())

e = 0.04
A=4*n/e*d**2/r_Z**2
vergleiche_aktivitaet(A, filename='3.2')

# Raumwinkelkorrektur
l_c = l/2.
A_k1=4*n/e*(d+l_c)**2/r_Z**2
k1 = A_k1/A
vergleiche_aktivitaet(A_k1, k=k1, index=1, filename='3.3')

l_c = l/4.
A_k2=4*n/e*(d+l_c)**2/r_Z**2
k2 = A_k2/A
vergleiche_aktivitaet(A_k2, k=k2, index=2, filename='3.4')

# Absorptionskorrektur
A_k3 = A_k2*unc.umath.exp(k*7.9*const.gram/const.centi**3*1.4*const.milli)
k3 = A_k3/A_k2
e_soll = 4*n/A_erw*k3*(d+l_c)**2/r_Z**2
vergleiche_aktivitaet(A_k3, e_soll=e_soll, index=3, filename='3.5')
print "e_soll:", (e_soll[0]+e_soll[1]+e_soll[2])/3.

plt.clf()
plt.title(u'Diagramm 3.3: Aktivitätsmessung des '+r'$\gamma$'+u'-Präparats in Abhängigkeit des Abstands')
xlim = np.array([0.04, 0.22])
plt.xlim(*xlim/const.centi)
plt.xlabel(u'Abstand Präparat-Zählrohr '+r'$d \, [cm]$')
plt.ylabel(u'Aktivität '+r'$A \, [\frac{Zerfaelle}{\mu s}]$')
papstats.plot_data(d/const.centi, A/const.mega, label=u'Originalmessung')
papstats.plot_data(d/const.centi, A_k1/const.mega, label=u'Raumwinkelkorrektur mit '+r'$\frac{l}{2}$')
papstats.plot_data(d/const.centi, A_k2/const.mega, label=u'Raumwinkelkorrektur mit '+r'$\frac{l}{4}$')
papstats.plot_data(d/const.centi, A_k3/const.mega, label=u'zusätzliche Absorptionskorrektur')
plt.axhline(A_erw.n/const.mega, color='black', lw=2, label='Erwartungswert $'+papstats.pformat(A_erw.n, prec=3, label='A_{erw}', unit='s^{-1}')+'$')
plt.legend()
papstats.savefig_a4('3.3.png')


#####
print "\n# 5 Absorption von a-Strahlung"
#####

p1, p2, N = np.loadtxt('4.txt', skiprows=1, unpack=True)

p = unp.uarray(p1+(p2-p1)/2.,(p2-p1)/2.)
n = unp.uarray(N, np.sqrt(N))/60.

def fit_const(p, n_0):
    return (p-p)+n_0

def fit_lin(p, n_0, b):
    return p*b+n_0

sl = [slice(0, 6), slice(7, 12), slice(15, None)]
fit = [fit_const, fit_lin, fit_const]
fit_result = []
for i in range(3):
    fit_result.append(papstats.curve_fit(fit[i], p[sl[i]], n[sl[i]]))

n_1, n_2 = fit_result[0][0][0], fit_result[2][0][0]
n_H = n_1-(n_1-n_2)/2.
p_H = (n_H-fit_result[1][0][0])/fit_result[1][0][1]
print n_H, p_H

p_0 = 1013.
s_0 = unc.ufloat(3.95, 0.05)*const.centi
s1 = p_H/p_0*s_0
print s1

s2 = 2.35/1.43*const.centi
s3 = 0.68*const.centi

s = s1+s2+s3
print s/const.centi

Ea = unc.ufloat(5.5, 0.1)
print "Energie der a-Strahlung:"
papstats.print_rdiff(Ea, unc.ufloat(5.48,0))

plt.clf()
plt.title('Diagramm 3.4: '+r'$\alpha$'+u'-Strahlung in Abhängigkeit des Luftdrucks im Abstand '+r'$s_0$')
ylim = [0, 240]
plt.ylim(*ylim)
plt.xlabel(u'Luftdruck '+r'$p \, [mbar]$')
plt.ylabel(u'Zählrate '+r'$n \, [\frac{Ereignisse}{s}]$')
papstats.plot_data(p, n, label='Messpunkte')
eq = ['n=n_A', 'n=n_0+b*p', 'n=n_E']
plabels = [['n_A'], ['n_0', 'b'], ['n_E']]
punits = [['s^{-1}'], ['s^{-1}', r'(s \times mbar)^{-1}'], ['s^{-1}']]
for i in range(3):
    papstats.plot_fit(fit[i], fit_result[i][0], fit_result[i][1], np.linspace(p[sl[i]][0].n, p[sl[i]][-1].n), eq=eq[i], plabels=plabels[i], punits=punits[i])
plt.axvline(p_H.n, color='black', lw=2, ls='dotted', label='Reichweite der '+r'$\alpha$'+u'-Strahlung:\n$'+papstats.pformat(n_H, label=r'n_{\frac{1}{2}}', unit='s^{-1}')+'$\n$'+papstats.pformat(p_H, label=r'p_{\frac{1}{2}}', unit='mbar')+'$')
plt.fill_betweenx(np.linspace(ylim[0], ylim[1]), p_H.n-p_H.s, p_H.n+p_H.s, color='green', alpha=0.2)
plt.legend()
papstats.savefig_a4('3.4.png')

