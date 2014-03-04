# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 255 - Röntgenspektrometer
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
import prettytable as pt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats


#####
# Konstanten
#####

d_LiF = 201.4*const.pico
Z_LiF = 3

h_erw = unc.ufloat(6.62606957e-34, 0)

#####
# Funktionen
#####

def y_bragg(b, n=1):
    global d
    return 2*d*unp.sin(b/360.*2.*const.pi)/n

def h_planck(y, U_B):
    return const.e*U_B*y/const.c


global d
d = d_LiF


#####
print "\n# a: Grenzwellenlänge und Plancksche Konstante aus LiF Spektrum"
#####

b, n = np.loadtxt('1.a.txt', skiprows=1, unpack=True)

n = unp.uarray(n, np.sqrt(n*5)/5)

# Untergrund
def fit_U(b, n_U):
    return b-b+n_U
sl_U = slice(0, 13)
popt_U, pstats_U = papstats.curve_fit(fit_U, b[sl_U], n[sl_U])
n_U = popt_U[0]
print "Untergrund:", papstats.pformat(n_U, format='.2u')

# Bremsspektrum-Fit mit Kramerscher Regel
def kramer(y, ymin, K):
    y = unp.nominal_values(y)
    return K*(y/ymin-1)/y**2
def fit_brems(b, ymin, K):
    return kramer(y=y_bragg(b), ymin=ymin, K=K)
sl_brems = ( (n <= 200) & ( (b <= 17) | (n <= 45) ) ) & (n > 20)
popt_brems, pstats_brems = papstats.curve_fit(fit_brems, b[sl_brems], n[sl_brems], p0=[4.133e-11, 1e-18])

# Extrapolation
def fit_lin(b, a, n_0):
    return a*b + n_0

sl_grenz=slice(13, 18)
popt_grenz, pstats_grenz = papstats.curve_fit(fit_lin, b[sl_grenz], n[sl_grenz])
b_G = (n_U-popt_grenz[1])/popt_grenz[0]
y_G = y_bragg(b_G)

# Berechnung
print "Grenzwellenlänge:", papstats.pformat(y_G/const.pico, label='y_G', unit='pm', format='.2u')
print "Plancksches Wirkungsquantum:"
papstats.print_rdiff(h_planck(y=y_G, U_B=30*const.kilo), h_erw)
# 2. Ordnung
print "2. Ordung ab:", papstats.pformat(unp.arcsin(y_G/d)/2/const.pi*360, label='b', unit='°', format='.2u')

# Plot
plt.clf()
plt.title(u'Diagramm 3.1: Bremsspektrum von LiF mit Untergrund und Extrapolation am kurzwelligen Ende')
axmain = plt.subplot(111)
plt.xlabel(ur'Bestrahlungswinkel $\beta \, [^\circ]$')
plt.ylabel(ur'Zählrate $n \, [\frac{Ereignisse}{s}]$')
xlim = [b[0], b[-1]]
xspace = np.linspace(*xlim, num=1000)
axmain.set_xlim(*xlim)
papstats.plot_data(b, n, label='Messpunkte')
plt.fill_between(unp.nominal_values(b), 0, unp.nominal_values(n), color='g', alpha=0.2)

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

zoomxlim = [3, 7]
zoomylim = [0, 200]
zoomxspace = np.linspace(*zoomxlim, num=100)
axzoom = zoomed_inset_axes(axmain, 2.5, loc=1)
axzoom.set_xlim(*zoomxlim)
axzoom.set_ylim(*zoomylim)
mark_inset(axmain, axzoom, loc1=4, loc2=2, fc="none", ec="0.5")
papstats.plot_data(b, n)
plt.fill_between(unp.nominal_values(b), 0, unp.nominal_values(n), color='g', alpha=0.2)
papstats.plot_fit(fit_U, popt_U, xspace=zoomxspace, eq=ur'n=n_U, \, \beta \in ['+str(b[sl_U][0])+','+str(b[sl_U][-1])+ur']^\circ', punits=['s^{-1}'], lw=2)
papstats.plot_fit(fit_lin, popt_grenz, xspace=zoomxspace, eq=r'n=a* \beta +n_0, \, \beta \in ['+str(b[sl_grenz][0])+','+str(b[sl_grenz][-1])+ur']^\circ', punits=[r'(^\circ s)^{-1}', 's^{-1}'], lw=2)

axmain.legend(loc='upper left')
axzoom.legend(loc='upper left')
papstats.savefig_a4('3.1.png')

    
#####
print "\n# b: Spektralanalyse mit LiF"
#####

def analyze_spektrallinien(fileprefix, figindex, crstl, sl, d=None, y=None):

    data = np.append(np.loadtxt(fileprefix+'.b.1.txt', skiprows=1), np.loadtxt(fileprefix+'.b.2.txt', skiprows=1), axis=0)

    b, n = data[:,0], data[:,1]
    n = unp.uarray(n, np.sqrt(n*20)/20)
    
    sl = [ [(b >= bounds[0]) & (b <= bounds[1]) for bounds in sl_row] for sl_row in sl]

    def fit_gauss(x, m, s, A, n_0):
        return A/np.sqrt(2*const.pi)/s*np.exp(-((x-m)**2)/2/(s**2))+n_0
    
    r = []
    
    plt.clf()
    papstats.plot_data(b,n)
    papstats.savefig_a4('3.'+str(figindex)+'.a.png')

    plt.clf()
    plt.suptitle('Diagramm 3.'+str(figindex)+u': Spektrallinien von Molybdän bei Vermessung mit einem '+crstl+'-Kristall')
    for i in range(2):
        r.append([])
        # Linie
        for k in range(2):
            # Ordnung
            b_k = b[sl[i][k]]
            n_k = n[sl[i][k]]
            xspace = np.linspace(b_k[0], b_k[-1], num=1000)
            plt.subplot(2,2,i*2+k+1)
            plt.xlim(xspace[0], xspace[-1])
            if i==1:
                plt.xlabel(u'Bestrahlungswinkel '+r'$\beta \, [^\circ]$')
            if k==0:
                plt.ylabel(u'Zählrate '+r'$n \, [\frac{Ereignisse}{s}]$')
            plt.title('$K_{'+(r'\alpha' if i==0 else r'\beta')+'}$ ('+str(k+1)+'. Ordnung)')
            papstats.plot_data(b_k, n_k)
            # Gauss-Fit
            popt, pstats = papstats.curve_fit(fit_gauss, b_k, n_k, p0=[b_k[0]+(b_k[-1]-b_k[0])/2, (b_k[-1]-b_k[0])/4, np.sum(n_k).n, n_k[0].n])
            plt.fill_between(b_k, 0, unp.nominal_values(n_k), color='g', alpha=0.2)
            FWHM = popt[1]*2*unp.sqrt(2*unp.log(2))
            plt.hlines(popt[3].n+(fit_gauss(xspace, *unp.nominal_values(popt)).max()-popt[3].n)/2, popt[0].n-FWHM.n/2, popt[0].n+FWHM.n/2, color='black', lw=2, label='$'+papstats.pformat(FWHM, label='FWHM', unit=r'^\circ')+'$')
            papstats.plot_fit(fit_gauss, popt, xspace=xspace, plabels=[r'\mu', r'\sigma', 'A', 'n_0'], punits=['^\circ', '^\circ', 's^{-1}', 's^{-1}'])
            plt.ylim(unp.nominal_values(n_k).min()-n_k[unp.nominal_values(n_k).argmin()].s, unp.nominal_values(n_k).max()+(unp.nominal_values(n_k).max()-unp.nominal_values(n_k).min()))
            plt.legend(loc='upper center', prop={'size':10})

            b_S = unc.ufloat(popt[0].n, np.abs(popt[1].n))
            print "Winkel:", papstats.pformat(b_S, unit='°', format='.2u')
            if y is None:
                r[i].append(y_bragg(b_S, n=k+1))
                print "Wellenlänge der Linie:", papstats.pformat(r[i][k]/const.pico, label='y', unit='pm', format='.2u')
            if d is None:
                r[i].append((k+1)*y[i][k]/unc.umath.sin(b_S*const.degree))
                print "Gitterkonstante:", papstats.pformat(r[i][k]/const.pico, label='a', unit='pm', format='.2u')

    papstats.savefig_a4('3.'+str(figindex)+'.png')

    return r

y = analyze_spektrallinien(fileprefix='1', figindex=2, crstl='LiF', sl=[[(9.4, 10.8), (19.7, 21.5)], [(8.4, 9.5), (17.4, 19)]], d=d_LiF)
y_m = np.mean(y, axis=1)
y_erw = np.array([71.1, 63.1])*const.pico
print "Wellenlängen der Linien (gemittelt):"
for i in range(len(y_m)):
    papstats.print_rdiff(y_m[i]/const.pico, y_erw[i]/const.pico)


#####
print "\n# c: Plancksches Wirkungsquantum aus Spannungsextrapolation"
#####

U, n = np.loadtxt('1.c.txt', skiprows=1, unpack=True)
U *= const.kilo

n = unp.uarray(n, np.sqrt(n*5)/5)
n -= n[0]

sl = slice(3, None)

def fit_lin(U, a, n_0):
    return a*U+n_0

popt, pstats = papstats.curve_fit(fit_lin, U[sl], n[sl])

U_E = -popt[1]/popt[0]
print "Einsetzspannung:", papstats.pformat(U_E, format='.2u')
d = 201.4*const.pico
b = 7.5/360.*2*const.pi
h = const.e*U_E*2*d*np.sin(b)/const.c
papstats.print_rdiff(h, unc.ufloat(6.626*1e-34, 0))
print (h-h_erw)*const.c/2*d*np.sin(b)/const.e

plt.clf()
plt.title('Diagramm 3.3: Extrapolation zur Bestimmung der Einsetzspannung')
plt.xlabel(u'Beschleunigungsspannung '+r'$U_B \, [kV]$')
plt.ylim(-50, 300)
papstats.plot_data(U/const.kilo, n)
papstats.plot_fit(fit_lin, popt, pstats, np.linspace(U[0], U[-1]), xscale=1./const.kilo, eq=r'n=a*U_B+n_0', punits=['(Vs)^{-1}', 's^{-1}'])
plt.axhline(0, color='black')
plt.legend(loc='upper left')
papstats.savefig_a4('3.3.png')


#####
print "\n# 2: Spektralanalyse mit NaCl"
#####

a = analyze_spektrallinien(fileprefix='2', figindex=4, crstl='NaCl', sl=[[(6.8, 7.7), (14, 15.3)], [(5.8, 6.9), (11, 13.5)]], y=y)
a = np.mean(a)
print "Gitterkonstante (gemittelt):", papstats.pformat(a/const.pico, unit='pm', format='.2u')

NA = 4*58.44*const.gram/(2.164*const.gram/const.centi**3)/a**3
NA_erw = 6.0221e23
print "Avogadrozahl:"
papstats.print_rdiff(NA/1e23, NA_erw/1e23)

