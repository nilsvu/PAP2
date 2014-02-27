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
import prettytable as pt

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats


#####
print "# 1 Untergrundbestimmung"
#####

data = np.loadtxt('1.dat')

N = data[:,1]

table = pt.PrettyTable()
table.add_column('t [s]', np.arange(len(N))*10+5, align='r')
table.add_column('N', N.astype(int), align='r')
with open("Resources/table_N_U.txt", "w") as text_file:
    text_file.write(table.get_string())

N_U0 = unc.ufloat(np.mean(N),np.std(N)/np.sqrt(len(N)))/10.
print N_U0


#####
# Zerfallsfunktionen
#####

def fit_silber(t, t_1, N_1, t_2, N_2):
    global N_U
    return N_U.n+N_1*np.exp(-t/t_1)+N_2*np.exp(-t/t_2)

def fit_indium(t, t_0, N_0):
    global N_U
    return N_U.n+N_0*np.exp(-t/t_0)


#####
print "\n# 2 Halbwertszeitbestimmung"
#####

def compute_hwz(N_list, ttor, fit, plotname, title, sl=slice(None,None), Uscale=1, p0=None, eq=None, plabels=None, punits=None, Th_erw=None):
    
    N = np.sum(unp.uarray(N_list,np.sqrt(N_list)), axis=0)
    t = np.arange(len(N))*ttor+ttor/2.

    table = pt.PrettyTable()
    table.add_column('t [s]', t.astype(int), align='r')
    if len(N_list) > 1:
        for i in range(len(N_list)):
            table.add_column('N'+str(i+1), N_list[i].astype(int), align='r')
        table.add_column('Summe', N, align='r')
    else:
        table.add_column('N', N, align='r')
    with open("Resources/table_"+plotname+".txt", "w") as text_file:
        text_file.write(table.get_string())


    global N_U
    N_U = N_U0*Uscale*ttor
    popt, pstats = papstats.curve_fit(fit, t[sl], N[sl], p0=p0)

    # Untergrundfehler
    N_U = (N_U0-N_U0.s)*Uscale*ttor
    popt_min, pstats_min = papstats.curve_fit(fit, t[sl], N[sl], p0=p0)
    N_U = (N_U0+N_U0.s)*Uscale*ttor
    popt_max, pstats_max = papstats.curve_fit(fit, t[sl], N[sl], p0=p0)
    N_U = N_U0*Uscale*ttor
    s_U = unp.nominal_values(((np.abs(popt-popt_min)+np.abs(popt-popt_max))/2.))
    s_corrected = np.sqrt(unp.std_devs(popt)**2 + s_U**2)
    popt_corrected = unp.uarray(unp.nominal_values(popt),s_corrected)
    
    # Halbwertszeit
    Th = popt_corrected[::2]*unc.umath.log(2)
    for i in range(len(Th)):
        papstats.print_rdiff(Th[i]/60, Th_erw[i]/60)

    # Plot
    plt.clf()
    plt.title('Diagramm '+plotname+': '+title)
    plt.xlabel('Messzeit $t \, [s]$')
    plt.ylabel('Ereigniszahl $N$')
    xspace = np.linspace(0, t[-1])
    papstats.plot_data(t, N, label='Messpunkte')
    papstats.plot_fit(fit, popt, pstats, xspace, eq=eq, plabels=plabels, punits=punits)
    plt.fill_between(xspace, fit(xspace, *unp.nominal_values(popt_min)), fit(xspace, *unp.nominal_values(popt_max)), color='g', alpha=0.2)
    Nmin = np.amin(unp.nominal_values(N))
    for i in range(len(Th)):
        plt.hlines(popt[1::2][i].n/2.+N_U.n, 0, Th[i].n, lw=2, label='Halbwertszeit $'+papstats.pformat(Th[i], label=r'T_{\frac{1}{2}}'+('^'+str(i+1) if len(Th)>1 else ''), unit='s')+'$')
    handles, labels = plt.gca().get_legend_handles_labels()
    p = plt.Rectangle((0, 0), 1, 1, color='g', alpha=0.2)
    handles.append(p)
    labels.append('Fit im '+r'$1 \sigma$'+'-Bereich von $N_U$:'+''.join(['\n$'+papstats.pformat(s_U[i], label='\Delta '+plabels[i]+'^{U}', unit=punits[i])+'$' for i in range(len(plabels))]))
    plt.legend(handles, labels)
    papstats.savefig_a4(plotname+'.png')


compute_hwz(N_list=[np.loadtxt('2.'+str(i+1)+'.dat')[:,1] for i in range(4)], ttor=10, Uscale=4, fit=fit_silber, plotname='3.1', title='Zerfall von Silber mit Untergrund', eq=r'N(t)=N_1*e^{-\frac{t}{\tau_1}}+N_2*e^{-\frac{t}{\tau_2}}+N_U', plabels=[r'\tau_1','N_1',r'\tau_2','N_2'], punits=['s',None,'s',None], Th_erw=unp.uarray([24.6, 144.6],0))

compute_hwz(N_list=[np.loadtxt('3.dat')[:,1]], sl=slice(1,None), ttor=120, fit=fit_indium, p0=[54*60/np.log(2),1000], plotname='3.2', title='Zerfall von Indium mit Untergrund', eq=r'N(t)=N_0*e^{-\frac{t}{\tau}}+N_U', plabels=[r'\tau','N_0'], punits=['s',None], Th_erw=unp.uarray([3240], 0))





