# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 256 - Röntgenfluoreszenz
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import scipy as sp
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
# Konstanten
#####

Ry = 13.6*const.eV

elements = [('Fe', 26), ('Mo', 42), ('Zr', 40), ('Zn', 30), ('Cu', 29), ('Ni', 28), ('Ag', 47), ('Ti', 22)]
legierungen = ['L1', 'L2', 'L3', 'L5']

#####
# Funktionen
#####
def moseley(Z, s=1, n1=1, n2=2):
    return (Z-s)*np.sqrt(n1**(-2)-n2**(-2))


#####
# Messdaten verarbeiten
#####
if False:
    data = np.append(np.loadtxt('Messung_Fe_Mo_Zr_Zn_Cu_Ni_ag.txt', skiprows=1), np.loadtxt('Messung_Fe_Mo_Ti_1_2_3_5.txt', skiprows=1, usecols=range(4, 14)), axis=1)
    filenames = list(np.array(elements, dtype=object)[:,0]) + legierungen
    for i in range(len(filenames)):
        np.savetxt('data/'+filenames[i]+'.txt', data[:,2*i:2*i+2], delimiter='\t', header='E [eV]\tn [1/s]')


#####
# Kombiniertes Spektrum
#####

data = np.array([np.loadtxt('data/'+e[0]+'.txt', skiprows=1) for e in elements], dtype=object)
data[:,:,1] = unp.uarray(data[:,:,1], data[:,:,1]/np.sqrt(180))

plt.clf()
plt.title(u'Diagramm 3.1: Röntgenfluoreszenzspektrum verschiedener Elemente')
plt.xlabel('Strahlungsenergie $E \, [keV]$')
plt.ylabel(u'Zählrate '+r'$n \, [\frac{Ereignisse}{s}]$')
for i in range(data.shape[0]):
    p = papstats.plot_data(data[i,:,0], data[i,:,1], label=elements[i][0]+' ('+str(elements[i][1])+')', elinewidth=0.5, capsize=4)
    plt.fill_between(data[i,:,0], 0, unp.nominal_values(data[i,:,1]), color=p[0].get_color(), alpha=0.1)
#plt.xlim(np.min(data[:,:,0]), np.max(data[:,:,0]))
#plt.ylim(np.min(unp.nominal_values(data[:,:,1])), np.max(unp.nominal_values(data[:,:,1])))
plt.xlim(0, 26)
plt.legend()
papstats.savefig_a4('3.1.png')


#####
# Spektrallinien
#####

Z = np.array([e[1] for e in elements])

converters = dict.fromkeys(range(1, 4), unc.ufloat_fromstr)
E = np.loadtxt('elements.txt', skiprows=1, usecols=range(1, 3), converters=converters, dtype=object)*const.kilo*const.eV

# Abschirmungskonstante
ERy = unp.sqrt(E/Ry)
n1 = 1
n2 = np.array([2., 3.])
s = Z[:, np.newaxis]-ERy/np.sqrt(n1**(-2)-n2**(-2))

# Table output
table = pt.PrettyTable()
table.add_column('Element', [e[0]+' ('+str(e[1])+')' for e in elements])
for k in range(np.shape(E)[1]):
    l = ['alpha', 'beta'][k]
    table.add_column('E_'+l+' [keV]', [('{:.2u}'.format(e) if e.s!=0 else '{:.2f}'.format(e.n)) for e in E[:,k]/const.kilo/const.eV])
    table.add_column('(E_'+l+'/Ry)^1/2', [('{:.2u}'.format(e) if e.s!=0 else '{:.2f}'.format(e.n)) for e in ERy[:,k]])
    table.add_column('sigma_'+l, [('{:.2u}'.format(e) if e.s!=0 else '{:.2f}'.format(e.n)) for e in s[:,k]])
with open('Resources/3.1.txt', 'w') as file:
    file.write(table.get_string())

# Plot
plt.clf()
plt.suptitle(u'Diagramm 3.2: Überprüfung des Moseley\'schen Gesetzes')
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
gs.update(hspace=0.1)
xspace = np.linspace(20, 50)
ax1 = plt.subplot(gs[0])
ax1.set_ylabel(r'$\sqrt{\frac{E}{R_Y}}$')
ax1.set_xlim(20,50)
ax2 = plt.subplot(gs[1], sharex=ax1)
plt.setp(ax1.get_xticklabels(), visible=False)
ax2.set_xlabel('Ordnungszahl $Z$')
ax2.set_ylabel(u'Abschirmungskonstante $\sigma$')
ax3 = ax1.twiny()
ax3.set_xlim(20,50)
ax3.set_xticks(np.sort(Z))
ax3.set_xticklabels(np.array(elements)[:,0][np.argsort(Z)])
for k in range(np.shape(E)[1]):
    l = 'K_'+[r'\alpha', r'\beta'][k]
    p = papstats.plot_data(Z, ERy[:,k], label='$'+l+'$ Messpunkte', ax=ax1)
    pstats = papstats.PAPStats(ydata=unp.nominal_values(ERy[:,k][unp.std_devs(ERy[:,k])>0]), ymodel=moseley(Z[unp.std_devs(ERy[:,k])>0], s=k+1, n1=1, n2=k+2), sigma=unp.std_devs(ERy[:,k])[unp.std_devs(ERy[:,k])>0], ddof=2)
    ax1.plot(xspace, moseley(xspace, s=k+1, n1=1, n2=k+2), ls='dotted', label='$'+l+'$ nach Moseley\'schem Gesetz\n'+pstats.legendstring(), color=p[0].get_color())
    papstats.plot_data(Z, s[:,k], ax=ax2, color=p[0].get_color())
    ax2.axhline(k+1, color=p[0].get_color())
    ax2.fill_between(np.sort(Z), k+1, unp.nominal_values(s[:,k])[np.argsort(Z)], alpha=0.15, color=p[0].get_color())
ax1.set_ylim(None, 45)
#ax2.set_ylim(0.3, 2.3)
ax1.legend(loc='lower right')
papstats.savefig_a4('3.2.png')


#####
# Legierungen
#####

converters = dict.fromkeys(range(1), unc.ufloat_fromstr)
data = np.array([np.loadtxt('legierungen/'+l+'.txt', converters=converters, dtype=object) for l in legierungen])*const.kilo*const.eV
table_strings = []
for i in range(data.shape[0]):
    E_L = data[i]
    table = pt.PrettyTable()
    table.add_column(u'Moegliche Legierung '+legierungen[i]+':', ['Abweichung [keV]', 'Sigma'])
    containing_elements = []
    for j in range(E.shape[0]):
        diffs = np.abs([E[j] - E_l for E_l in E_L])/const.kilo/const.eV
        argmin = np.argmin(diffs)
        diff_min = diffs.flatten()[argmin]
        sigma = diff_min.n/diff_min.s
        if sigma <=1:
            containing_elements.append((j, sigma))
            table.add_column(elements[j][0],  np.array(papstats.rdiff(diff_min))[[0,2]])
    table_strings.append(table.get_string())
with open('Resources/3.2.txt', 'w') as file:
    file.write('\n'.join(table_strings))
