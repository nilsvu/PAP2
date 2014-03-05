# -*- coding: utf-8 -*-

'''
    Auswertung des Versuchs 255 - Röntgenfluoreszenz
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

elements = [('Fe', 26), ('Mo', 42), ('Zr', 40), ('Zn', 30), ('Cu', 29), ('Ni', 28), ('Ti', 22)]
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
# Messdaten anzeigen
#####

plt.clf()
data = np.array([np.loadtxt('data/'+e[0]+'.txt', skiprows=1) for e in elements])
for i in range(data.shape[0]):
    p = papstats.plot_data(data[i,:,0], data[i,:,1], label=elements[i][0])
    plt.fill_between(data[i,:,0], 0, data[i,:,1], color=p[0].get_color(), alpha=0.1)
plt.legend()
papstats.savefig_a4('3.1.png')


#####
# Spektrallinien
#####

converters = dict.fromkeys(range(1, 4), unc.ufloat_fromstr)
data = np.loadtxt('elements.txt', skiprows=1, usecols=range(1, 3), converters=converters, dtype=object)
elements = np.genfromtxt('elements.txt', dtype='str', usecols=range(1), skiprows=1)

Z = np.array(data[:,0], dtype=int)
E = data[:,1:]*const.kilo*const.eV
ERy = unp.sqrt(E/Ry)

n1 = 1
n2 = np.array([2., 3.])

s = np.transpose(Z-np.transpose(ERy/np.sqrt(n1**(-2)-n2**(-2))))

table = pt.PrettyTable()
table.add_column('Elemente', elements)
for k in range(np.shape(E)[1]):
    table.add_column('E'+str(k+1), E[:,k])
    table.add_column('(E'+str(k+1)+'/Ry)^1/2', ERy[:,k])
    table.add_column('sigma'+str(k+1), s[:,k])
print table

plt.clf()
xspace = np.linspace(20, 50)
for k in range(np.shape(E)[1]):
    papstats.plot_data(Z, unp.nominal_values(ERy[:,k]))
    plt.plot(xspace, moseley(xspace, s=k+1, n1=1, n2=k+2), ls='dotted')
papstats.savefig_a4('3.2.png')

plt.clf()
for k in range(np.shape(E)[1]):
    papstats.plot_data(Z, unp.nominal_values(s[:,k]))
    plt.axhline(k+1)
plt.ylim(0, 3)
papstats.savefig_a4('3.2.png')


converters = dict.fromkeys(range(3), unc.ufloat_fromstr)
print np.loadtxt('legierungen.txt', converters=converters)
