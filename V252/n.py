# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 252 - Aktivierung von Indium und von Silber mit thermischen Neutronen
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as const

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats


#####
print "# 1 Untergrundbestimmung"
#####

data = np.loadtxt('1.dat')

N = data[:,1]
N_U = unc.ufloat(np.mean(N),np.std(N))


# Zerfallsfunktion
def fit_silber(t, t_1, N_1, t_2, N_2):
    return N_U.n+N_1*np.exp(-t/t_1)+N_2*np.exp(-t/t_2)

def fit_indium(t, t_0, N_0):
    return N_U.n+N_0*np.exp(-t/t_0)


#####
print "\n# 2 Halbwertszeitbestimmung"
#####

def compute_hwz(N, ttor, fit, plotname, title, sl=slice(None,None), p0=None, eq=None, plabels=None, punits=None):
    
    N = unp.uarray(N,np.sqrt(N))
    t = np.arange(len(N))*ttor+ttor/2.

    popt, pstats = papstats.curve_fit(fit, t[sl], N[sl], p0=p0)

    plt.clf()
    plt.title('Diagramm '+plotname+': '+title)
    plt.xlabel('Messzeit $t \, [s]$')
    plt.ylabel('Ereigniszahl $N$')
    xspace = np.linspace(0, t[-1])
    papstats.plot_data(t,N)
    papstats.plot_fit(fit, popt, pstats, xspace, eq=eq, plabels=plabels, punits=punits)
    plt.legend()
    papstats.savefig_a4(plotname+'.png')


compute_hwz(N=np.sum([np.loadtxt('2.'+str(i+1)+'.dat')[:,1] for i in range(4)], axis=0), ttor=10, fit=fit_silber, plotname='3.1', title='Zerfall von Silber mit Untergrund', eq=r'N(t)=N_1*e^{-\frac{t}{\tau_1}}+N_2*e^{-\frac{t}{\tau_2}}+N_U', plabels=[r'\tau_1','N_1',r'\tau_2','N_2'], punits=['s',None,'s',None])

compute_hwz(N=np.loadtxt('3.dat')[:,1], sl=slice(1,None), ttor=120, fit=fit_indium, p0=[54*60/np.log(2),1000], plotname='3.2', title='Zerfall von Indium mit Untergrund', eq=r'N(t)=N_0*e^{-\frac{t}{\tau}}+N_U', plabels=[r'\tau','N_0'], punits=['s',None])





