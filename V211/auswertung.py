# -*- coding: utf-8 -*-

'''
Auswertung des Versuches 211 - Gekoppelte Pendel
Phyiskalisches Anfänger Praktikum II - Universität Heidelberg

Autoren: Christian Kohlstedde und Nils Fischer
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

# Unicode Hack
import sys
import codecs
#sys.stdout=codecs.getwriter('utf-8')(sys.stdout)


# Functions

def print_w(T1, T2):
	print papstats.pformat(T1, unit='s', label='T_1')
	print papstats.pformat(T2, unit='s', label='T_2')
	w1 = 2*const.pi/T1
	w2 = 2*const.pi/T2
	w = np.mean([w1, w2])
	print papstats.pformat(w1, unit='Hz', label='w_1')
	print papstats.pformat(w2, unit='Hz', label='w_2')
	print papstats.pformat(w, unit='Hz', label='w')
	return (w1, w2, w)


def schwingung_table(l, t1, t2, n):
	T1, T2 = t1 / n, t2 / n
	T = np.mean([T1, T2], axis=0)
	w = 2 * const.pi / T
	return (papstats.table(labels=['Kopplung', 'T_1', 'T_2', 'T', 'w'], units=[None, 's', 's', 's', 'Hz'], columns=[papstats.pformat(l / const.centi, label='l', unit='cm', format='c'), T1, T2, T, w]), w)


def vergleich_table(l, w, w_erw):
	diff = w - w_erw
	return papstats.table(labels=['Kopplung', 'w', 'w_erw', 'w-w_erw', 'Sigmabereich'], units=[None, 'Hz', 'Hz', 'Hz', None], columns=[papstats.pformat(l / const.centi, label='l', unit='cm', format='c'), w, w_erw, diff, unp.nominal_values(np.abs(diff))/unp.std_devs(diff)])


# Konstanten

dt = 0.05


print '\n# 2.1 Offsets'

data = np.loadtxt('offset.txt', skiprows=1)
offset1 = np.mean(data[:,1])
offset2 = np.mean(data[:,2])
print papstats.pformat(offset1, label='Offset phi_1', unit='grad')
print papstats.pformat(offset2, label='Offset phi_2', unit='grad')


print "\n# 2.2 Eigenfrequenzen der ungekoppelten Pendel"

t1 = unc.ufloat(25.51, dt) - unc.ufloat(1.27, dt)
t2 = unc.ufloat(25.17, dt) - unc.ufloat(0.88, dt)
n = 15
T1, T2 = t1 / n, t2 / n
print_w(T1, T2)


print "\n# 2.3 Gekoppelte Pendel"

data_symm = []
data_anti = []

l = unc.ufloat(19.2, 0.1) * const.centi
t1 = unc.ufloat(33.08, dt) - unc.ufloat(0.77, dt)
t2 = unc.ufloat(33.06, dt) - unc.ufloat(0.82, dt)
n = 20
data_symm.append([l, t1, t2, n])

t1 = unc.ufloat(32.19, dt) - unc.ufloat(1.17, dt)
t2 = unc.ufloat(31.42, dt) - unc.ufloat(0.44, dt)
n = 20
data_anti.append([l, t1, t2, n])

l = unc.ufloat(0.292, 0.001)
t1 = unc.ufloat(32.50, dt) - unc.ufloat(0.24, dt)
t2 = unc.ufloat(32.48, dt) - unc.ufloat(0.23, dt)
n = 20
data_symm.append([l, t1, t2, n])

t1 = unc.ufloat(33.10, dt) - unc.ufloat(0.77, dt)
t2 = unc.ufloat(33.84, dt) - unc.ufloat(1.48, dt)
n = 22
data_anti.append([l, t1, t2, n])

l = unc.ufloat(0.392, 0.001)
t1 = unc.ufloat(33.88, dt) - unc.ufloat(1.67, dt)
t2 = unc.ufloat(33.89, dt) - unc.ufloat(1.65, dt)
n = 20
data_symm.append([l, t1, t2, n])

t1 = unc.ufloat(28.36, dt) - unc.ufloat(1.05, dt)
t2 = unc.ufloat(29.03, dt) - unc.ufloat(1.73, dt)
n = 20
data_anti.append([l, t1, t2, n])

print "Symmetrische Schwingung:"
(table, w_1) = schwingung_table(*np.transpose(data_symm))
print table
print "Antisymmetrische Schwingung:"
(table, w_2) = schwingung_table(*np.transpose(data_anti))
print table


print "\n# 2.4 Superposition der Anregung"

data_schwing = []
data_schweb = []

l = unc.ufloat(0.392, 0.001)
t1 = unc.ufloat(29.30, dt) - unc.ufloat(1.87, dt)
t2 = unc.ufloat(34.10, dt) - unc.ufloat(6.68, dt)
n = 20
data_schwing.append([l, t1, t2, n])

dt_Schwebung = 0.2
t1 = unc.ufloat(29.30, dt_Schwebung) - unc.ufloat(1.87, dt_Schwebung)
t2 = unc.ufloat(34.10, dt_Schwebung) - unc.ufloat(6.68, dt_Schwebung)
n = 1.5
data_schweb.append([l, t1, t2, n])

l = unc.ufloat(0.292, 0.001)
t1 = unc.ufloat(39.11, dt) - unc.ufloat(2.20, dt)
t2 = unc.ufloat(50.27, dt) - unc.ufloat(13.34, dt)
n = 25
data_schwing.append([l, t1, t2, n])

dt_Schwebung = 0.5
t1 = unc.ufloat(64.92, dt_Schwebung) - unc.ufloat(14.06, dt_Schwebung)
t2 = unc.ufloat(56.25, dt_Schwebung) - unc.ufloat(5.39, dt_Schwebung)
n = 1.5
data_schweb.append([l, t1, t2, n])

l = unc.ufloat(0.192, 0.001)
t1 = unc.ufloat(26.00, dt) - unc.ufloat(2.29, dt)
t2 = unc.ufloat(45.34, dt) - unc.ufloat(21.61, dt)
n = 15
data_schwing.append([l, t1, t2, n])

dt_Schwebung = 0.5
t1 = unc.ufloat(76.52, dt_Schwebung) - unc.ufloat(36.71, dt_Schwebung)
t2 = unc.ufloat(56.63, dt_Schwebung) - unc.ufloat(17.00, dt_Schwebung)
n = 0.5
data_schweb.append([l, t1, t2, n])

data_schwing = np.array(data_schwing)
data_schwing = data_schwing[np.argsort(data_schwing[:,0])]
data_schweb = np.array(data_schweb)
data_schweb = data_schweb[np.argsort(data_schweb[:,0])]

l = data_schwing[:,0]

print "Schwingung:"
(table, w_I) = schwingung_table(*np.transpose(data_schwing))
print table
w_I_erw = (w_1 + w_2) / 2.
print vergleich_table(l, w_I, w_I_erw)
print "Schwebung:"
(table, w_II) = schwingung_table(*np.transpose(data_schweb))
print table
w_II_erw = (w_2 - w_1) / 2.
print vergleich_table(l, w_II, w_II_erw)


print u"\n# Kopplungsgrade"

V_l = l**2 / np.roll(l, 1)**2
K = (w_2**2-w_1**2)/2/w_1**2
print K
V_w = K / np.roll(K, 1)
diff = V_l-V_w
K_labels = ['K_1', 'K_2', 'K_3']
K_labels = [K_labels[i] + '/' + np.roll(K_labels, 1)[i] for i in range(len(K_labels))]
print papstats.table(labels=['', 'V_l', 'V_w', 'V_l-V_w', 'Sigmabereich'], columns=[K_labels, V_l, V_w, diff, unp.nominal_values(np.abs(diff))/unp.std_devs(diff)])



