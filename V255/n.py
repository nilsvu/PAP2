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


#####
# Funktionen
#####

def y_bragg(b, n=1):
    global d
    return 2*d*unp.sin(b/360.*2.*const.pi)/n

def h_planck(y, U_B):
    return const.e*U_B*y/const.c



def analyze_spektrum(fileprefix, sl_U, sl_grenz):

    #####
    print "\n# a: Grenzwellenlänge und Plancksche Konstante aus Spektrum"
    #####

    b, n = np.loadtxt(fileprefix+'.a.txt', skiprows=1, unpack=True)

    n = unp.uarray(n, np.sqrt(n*5)/5)

    # Untergrund
    def fit_U(b, n_U):
        return b-b+n_U
    popt_U, pstats_U = papstats.curve_fit(fit_U, b[sl_U], n[sl_U])
    n_U = popt_U[0]
    print "Untergrund:", papstats.pformat(n_U, format='.2u')

    # Bremsspektrum-Fit mit Kramerscher Regel
    def kramer(y, ymin, K):
        return K*(y/ymin-1)/y**2
    def fit_brems(b, ymin, K, b_0):
        global d
        return kramer(y=y_bragg(b, d), ymin=ymin, K=K)+b_0

    #sl_brems = ( (n <= 200) & ( (b <= 17) | (n <= 45) ) ) & (n > 20)

    #popt_brems, pstats_brems = papstats.curve_fit(fit_brems, b[sl_brems], n[sl_brems], p0=[4.133e-11, 1e-18, 0])

    # Extrapolation
    def fit_lin(b, a, n_0):
        return a*b + n_0

    popt_grenz, pstats_grenz = papstats.curve_fit(fit_lin, b[sl_grenz], n[sl_grenz])

    y_G = y_bragg((n_U-popt_grenz[1])/popt_grenz[0])


    # Berechnung

    print "Grenzwellenlänge:", papstats.pformat(y_G/const.pico, label='y_G', unit='pm', format='.2u')
    print "Plancksches Wirkungsquantum:", papstats.pformat(h_planck(y=y_G, U_B=30*const.kilo), label='h', unit='Js', format='.2u')
    # 2. Ordnung
    print "2. Ordung ab:", papstats.pformat(unp.arcsin(y_G/d)/2/const.pi*360, label='b', unit='°', format='.2u')


    # Plot
    plt.clf()
    plt.title(u'Diagramm 3.'+fileprefix+'.1')
    axmain = plt.subplot(111)
    plt.xlabel(ur'Bestrahlungswinkel $\beta \, [°]$')
    plt.ylabel(ur'Zählrate $n \, [\frac{Ereignisse}{s}]$')
    xlim = [-2, 23]
    #ylim = [0, 3e2]
    xspace = np.linspace(*xlim, num=1000)
    axmain.set_xlim(*xlim)
    #axmain.set_ylim(*ylim)
    papstats.plot_data(b, n, label='Messpunkte')
    #plt.plot(unp.nominal_values(b[sl_brems]), unp.nominal_values(n[sl_brems]))
    #plt.plot(xspace, fit_brems(xspace, 4.133e-11, 1e-18, 0))
    papstats.plot_fit(fit_U, popt_U, pstats_U, xspace)
    #papstats.plot_fit(fit_brems, popt_brems, pstats_brems, xspace)

    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
    from mpl_toolkits.axes_grid1.inset_locator import mark_inset

    zoomxlim = [5, 7]
    zoomylim = [0, 200]
    zoomxspace = np.linspace(*zoomxlim, num=100)
    axzoom = zoomed_inset_axes(axmain, 3, loc=2)
    axzoom.set_xlim(*zoomxlim)
    axzoom.set_ylim(*zoomylim)
    #plt.xticks(visible=False)
    plt.yticks(visible=False)
    mark_inset(axmain, axzoom, loc1=1, loc2=3, fc="none", ec="0.5")
    papstats.plot_data(b, n)
    papstats.plot_fit(fit_U, popt_U, pstats_U, zoomxspace)
    papstats.plot_fit(fit_lin, popt_grenz, pstats_grenz, zoomxspace)

    axmain.legend()
    papstats.savefig_a4('3.'+fileprefix+'.1.png')

    
def analyze_spektrallinien(fileprefix, sl, d=None, y=None):

    #####
    print "\n# b: Spektrallinienanalyse"
    #####

    data = np.append(np.loadtxt(fileprefix+'.b.1.txt', skiprows=1), np.loadtxt(fileprefix+'.b.2.txt', skiprows=1), axis=0)

    b, n = data[:,0], data[:,1]
    n = unp.uarray(n, np.sqrt(n*20)/20)
    
    sl = [ [(b > bounds[0]) & (b < bounds[1]) for bounds in sl_row] for sl_row in sl]

    def fit_gauss(x, m, s, A, n_0):
        return A/np.sqrt((s**2)*2*const.pi)*np.exp(-((x-m)**2)/2/(s**2))+n_0
    
    r = []

    plt.clf()
    for i in range(2):
        r.append([])
        # Linie
        for k in range(2):
            # Ordnung
            b_k = b[sl[i][k]]
            n_k = n[sl[i][k]]
            xspace = np.linspace(b_k[0], b_k[-1], num=1000)
            plt.subplot(1,4,i*2+k+1)
            plt.title('$K_{'+(r'\alpha' if i==0 else r'\beta')+'}$ ('+str(k+1)+'. Ordnung)')
            papstats.plot_data(b_k, n_k)
            
            # Gauss-Fit
            popt, pstats = papstats.curve_fit(fit_gauss, b_k, n_k, p0=[b_k[0]+(b_k[-1]-b_k[0])/2, (b_k[-1]-b_k[0])/3, 1, n_k[0].n])
            papstats.plot_fit(fit_gauss, popt, pstats, xspace)

            if y is None:
                r[i].append(y_bragg(popt[0], n=k+1))
                print "Wellenlänge der Linie:", papstats.pformat(r[i][k]/const.pico, label='y', unit='pm', format='.2u')
            if d is None:
                r[i].append((k+1)*y[i]/unc.umath.sin(popt[0]/360*2*const.pi))
                print "Gitterkonstante:", papstats.pformat(r[i][k], label='a', unit='m', format='.2u')

    papstats.savefig_a4('3.'+fileprefix+'.2.png')

    return r


#####
print "### LiF"
#####

d = d_LiF
analyze_spektrum(fileprefix='1', sl_U=slice(0, 13), sl_grenz=slice(13, 17))
y = analyze_spektrallinien(fileprefix='1', sl=[[(8.3, 9.5), (17.25, 18.75)], [(9.5, 10.6), (19.8, 21.5)]], d=d_LiF)
y = np.mean(y, axis=0)
print "Wellenlängen der Linien (gemittelt):", y


#####
print "### NaCl"
#####

d = analyze_spektrallinien(fileprefix='2', sl=[[(6, 6.8), (12, 14)], [(6.8, 7.8), (14, 15)]], y=y)
d = np.mean(d)
print "Gitterkonstante (gemittelt):", d

NA = 58.44*const.gram/2/(2.164*const.gram/const.centi)/d**3
print "Avogadrozahl:", NA


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

U_0 = -popt[1]/popt[0]
d = 201.4*const.pico
b = 7.5/360.*2*const.pi
h = const.e*U_0*2*d*np.sin(b)/const.c
papstats.print_rdiff(h, unc.ufloat(6.626*1e-34, 0))

plt.clf()
papstats.plot_data(U, n)
papstats.plot_fit(fit_lin, popt, pstats, np.linspace(U[0], U[-1]))
plt.legend()
papstats.savefig_a4('3.3.png')
