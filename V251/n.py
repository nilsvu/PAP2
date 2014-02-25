# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 251 - Statistik des radioaktiven Zerfalls
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.special import gamma

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats

if False:
    ydata, xdata = np.histogram(np.random.normal(10,2,size=10000), bins=1000)
    xdata = xdata[:-1]+np.diff(xdata)/2
    xdata = xdata[ydata>0]
    ydata = unp.uarray(ydata, np.sqrt(ydata))
    ydata = ydata[ydata>0]

    def fit_gauss(x, m, s, A):
        return A/np.sqrt((s**2)*2*const.pi)*np.exp(-((x-m)**2)/2./(s**2))

    popt, pstats = papstats.curve_fit(fit_gauss, xdata, ydata)

    plt.clf()
    papstats.plot_data(xdata, ydata)
    papstats.plot_fit(fit_gauss, popt, pstats, np.linspace(0,20))
    plt.legend()
    papstats.savefig_a4('0.png')

    import sys
    sys.exit()



#####
print('# 1 (Plateaubereich des Zählrohrs)')
#####

data = np.loadtxt('2.txt', skiprows=1)

U = unp.uarray(data[:,0], 10)
N = data[:,1]
N = unp.uarray(N, np.sqrt(N))/30

def fit_platlin(x, c, N_0):
    return c*x+N_0
def fit_platconst(x, N_0):
    return np.zeros(len(x))+N_0

popt_const, pstats_const = papstats.curve_fit(fit_platconst, U[1:], N[1:], p0=[60])
popt_lin, pstats_lin = papstats.curve_fit(fit_platlin, U[1:], N[1:])
popt_lin2, pstats_lin2 = papstats.curve_fit(fit_platlin, U[8:], N[8:])

plt.clf()
plt.title(u'Diagramm 3.1: Vermessung des Plateaubereichs der Zählrohrkennlinie')
papstats.plot_data(U, N)
papstats.plot_fit(fit_platconst, popt_const, pstats_const, np.linspace(U[1].n, U[-1].n), eq='N=N_0', ls='dashed', lw=2)
papstats.plot_fit(fit_platlin, popt_lin, pstats_lin, np.linspace(U[1].n, U[-1].n), eq='N=c*U_Z+N_0', ls='dotted', lw=2)
papstats.plot_fit(fit_platlin, popt_lin2, pstats_lin2, np.linspace(U[8].n, U[-1].n), eq='N=c*U_Z+N_0, \, U_Z \in [600,700]V', lw=2)
plt.xlabel(u'Zählrohrspannung '+r'$U_Z \, [V]$')
plt.ylabel(u'Zählrate '+r'$\frac{N}{t} \, [\frac{Ereignisse}{s}]$')
plt.xlim(430,720)
plt.ylim(10,65)
plt.legend(loc='lower right')
papstats.savefig_a4('3.1.png')

#####
print('\n# 2 (Untersuchung des Plateauanstiegs)')
#####

data = np.loadtxt('3.txt', skiprows=1)

t = data[:,0]
N1 = data[:,1]
N1 = unp.uarray(N1,np.sqrt(N1))
N2 = data[:,2]
N2 = unp.uarray(N2,np.sqrt(N2))

Ndiff = N1-N2
papstats.print_rdiff(N1[0],N2[0])
papstats.print_rdiff(N1[1],N2[1])

# Zählrate
Z1 = N1/t
Z2 = N2/t
T = (Z1+Z2)/(Z1-Z2)**2*1e4
print T/(60*24)


#####
print('\n# 4 (Verifizierung der statistischen Natur des radioaktiven Zerfalls)')
#####

# Gaussverteilung

def fit_gauss(x, m, s, A):
    return A/np.sqrt((s**2)*2*const.pi)*np.exp(-((x-m)**2)/2/(s**2))

# Poissonverteilung

def fit_poisson(x, m, A):
    return A*np.exp(-m)*(m**x)/gamma(x+1)


def compare_gauss_poisson(t, data, p0, title, filename, xlim, ylim):

    N = data[:,0]
    n = data[:,1]
    n = unp.uarray(n,np.sqrt(n))

    sl = (n >= 10) # TODO: Häufigkeit n mindestens 10

    # Fit

    popt_gauss, pstats_gauss = papstats.curve_fit(fit_gauss, N[sl], n[sl], p0=p0, sigma=unp.std_devs(n[sl]))

    popt_poisson, pstats_poisson = papstats.curve_fit(fit_poisson, N[sl], n[sl], p0=[p0[0],p0[2]], sigma=unp.std_devs(n[sl]))

    # Plot
    
    for log in [False, True]:
        plt.clf()
        plt.title('Diagramm '+filename+('.b' if log else '.a')+': '+title + (' (logarithmisch)' if log else ''))
        if log:
            plt.yscale('log')
        papstats.plot_data(N/t, n)
        xrange = 4*popt_gauss[1].n
        xspace = np.linspace(xlim[2 if log else 0]*t,xlim[3 if log else 1]*t,num=200)
        papstats.plot_fit(fit_gauss, popt_gauss, pstats_gauss, xspace, xscale=1./t, eq=r'G(N;\mu,\sigma)', plabels=[r'\mu',r'\sigma','A'])
        papstats.plot_fit(fit_poisson, popt_poisson, pstats_poisson, xspace, xscale=1./t, eq=r'P(N;\mu)', plabels=[r'\mu','A'], ls='dashed')
        plt.xlim(xspace[0]/t,xspace[-1]/t)
        plt.ylim(ylim[2 if log else 0],ylim[3 if log else 1])
        plt.xlabel(u'Zählrate '+r'$Z=\frac{N}{t} \, [\frac{Ereignisse}{s}]$')
        plt.ylabel(u'Häufigkeit '+r'$n$')
        plt.legend(loc=('lower center' if log else 'upper right'))
        papstats.savefig_a4(filename+('.b' if log else '.a')+'.png')
    
    # Residuum
    plt.clf()
    plt.title('Diagramm '+filename+'.c: Residuum')
    plt.hist(fit_gauss(unp.nominal_values(N), *unp.nominal_values(popt_gauss))-unp.nominal_values(n), bins=30)
    plt.hist(pstats_gauss.residual, bins=30)
    plt.hist(pstats_poisson.residual, bins=30)
    papstats.savefig_a4(filename+'.c.png')


compare_gauss_poisson(t=0.5, data=np.loadtxt('4.dat'), p0=[57.768,7.659,2094], title=u'Vergleich der Poisson- und Gauß-Verteilung bei großen Zählraten', filename='3.2', xlim=[60,200,60,180], ylim=[-10,130,1e-1,1.5e2])

compare_gauss_poisson(t=0.1, data=np.loadtxt('5.dat'), p0=[4.134,2.050,5246], title=u'Vergleich der Poisson- und Gauß-Verteilung bei kleinen Zählraten', filename='3.3', xlim=[-10,140,-5,140], ylim=[-100,1100,1e-2,1.5e3])





