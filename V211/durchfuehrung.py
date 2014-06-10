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


print '### 2.1 Offsets'
data = np.loadtxt('1.offset.txt', skiprows=3)

offset1 = np.mean(data[:,1])
offset2 = np.mean(data[:,2])

print u'Offset des Pendel 1', offset1
print u'Offset des Pendel 2', offset2


print "### 2.2 Eigenfrequenzen der ungekoppelten Pendel ###"

t1 = unc.ufloat(25.51, 0.05) - unc.ufloat(1.27, 0.05)
t2 = unc.ufloat(25.17, 0.05) - unc.ufloat(0.88, 0.05)

n = 15

T1, T2 = t1 / n, t2 / n

print papstats.pformat(t1, unit='s', label='t_1')
print papstats.pformat(t2, unit='s', label='t_2')
print papstats.pformat(T1, unit='s', label='T_1')
print papstats.pformat(T2, unit='s', label='T_2')

print ''
print "### 2.3 Gekoppelte Pendel ###"

l = unc.ufloat(0.192, 0.001)
print papstats.pformat(l, unit='m', label=u'Länge l')

print "symmetrische Eigenschwingungen"
t1 = unc.ufloat(33.08, 0.05) - unc.ufloat(0.77, 0.05)
t2 = unc.ufloat(33.06, 0.05) - unc.ufloat(0.82, 0.05)

n = 20

T1, T2 = t1 / n, t2 / n

print papstats.pformat(t1, unit='s', label='t_1')
print papstats.pformat(t2, unit='s', label='t_2')
print papstats.pformat(T1, unit='s', label='T_1')
print papstats.pformat(T2, unit='s', label='T_2')


print "antisymmetrische Eigenschwingungen"
t1 = unc.ufloat(32.19, 0.05) - unc.ufloat(1.17, 0.05)
t2 = unc.ufloat(31.42, 0.05) - unc.ufloat(0.44, 0.05)

n = 20

T1, T2 = t1 / n, t2 / n

print papstats.pformat(t1, unit='s', label='t_1')
print papstats.pformat(t2, unit='s', label='t_2')
print papstats.pformat(T1, unit='s', label='T_1')
print papstats.pformat(T2, unit='s', label='T_2')

print ''

l = unc.ufloat(0.292, 0.001)
print papstats.pformat(l, unit='m', label=u'Länge l')


print "symmetrische Eigenschwingungen"
t1 = unc.ufloat(32.50, 0.05) - unc.ufloat(0.24, 0.05)
t2 = unc.ufloat(32.48, 0.05) - unc.ufloat(0.23, 0.05)

n = 20

T1, T2 = t1 / n, t2 / n

print papstats.pformat(t1, unit='s', label='t_1')
print papstats.pformat(t2, unit='s', label='t_2')
print papstats.pformat(T1, unit='s', label='T_1')
print papstats.pformat(T2, unit='s', label='T_2')

print "antisymmetrische Eigenschwingungen"
t1 = unc.ufloat(33.10, 0.05) - unc.ufloat(0.77, 0.05)
t2 = unc.ufloat(33.84, 0.05) - unc.ufloat(1.48, 0.05)

n = 22

T1, T2 = t1 / n, t2 / n

print papstats.pformat(t1, unit='s', label='t_1')
print papstats.pformat(t2, unit='s', label='t_2')
print papstats.pformat(T1, unit='s', label='T_1')
print papstats.pformat(T2, unit='s', label='T_2')

print ''

l = unc.ufloat(0.392, 0.001)
print papstats.pformat(l, unit='m', label=u'Länge l')

print "symmetrische Eigenschwingungen"
t1 = unc.ufloat(33.88, 0.05) - unc.ufloat(1.67, 0.05)
t2 = unc.ufloat(33.89, 0.05) - unc.ufloat(1.65, 0.05)

n = 20

T1, T2 = t1 / n, t2 / n

print papstats.pformat(t1, unit='s', label='t_1')
print papstats.pformat(t2, unit='s', label='t_2')
print papstats.pformat(T1, unit='s', label='T_1')
print papstats.pformat(T2, unit='s', label='T_2')

print "antisymmetrische Eigenschwingungen"
t1 = unc.ufloat(28.36, 0.05) - unc.ufloat(1.05, 0.05)
t2 = unc.ufloat(29.03, 0.05) - unc.ufloat(1.73, 0.05)

n = 20

T1, T2 = t1 / n, t2 / n

print papstats.pformat(t1, unit='s', label='t_1')
print papstats.pformat(t2, unit='s', label='t_2')
print papstats.pformat(T1, unit='s', label='T_1')
print papstats.pformat(T2, unit='s', label='T_2')

print ''
print "### 2.4 Superposition der Anrengung ###"

print ''
print "## Messung 1"

l = unc.ufloat(0.392, 0.001)
print u"Länge", papstats.pformat(l)

print "Schwingungsperiodendauer"
dt_Schwingung = 0.05
t1_1 = unc.ufloat(29.30, dt_Schwingung) - unc.ufloat(1.87, dt_Schwingung)
t2_1 = unc.ufloat(34.10, dt_Schwingung) - unc.ufloat(6.68, dt_Schwingung)

n = 20

T1_1, T2_1 = t1_1 / n, t2_1 / n

print papstats.pformat(t1_1, unit='s', label='t1_1')
print papstats.pformat(t2_1, unit='s', label='t2_1')
print papstats.pformat(T1_1, unit='s', label='T1_1')
print papstats.pformat(T2_1, unit='s', label='T2_1')

print "Schwebungsperiodendauer"
dt_Schwebung = 0.2
t1_2 = unc.ufloat(29.30, dt_Schwebung) - unc.ufloat(1.87, dt_Schwebung)
t2_2 = unc.ufloat(34.10, dt_Schwebung) - unc.ufloat(6.68, dt_Schwebung)

n = 1.5

T1_2, T2_2 = t1_2 / n, t2_2 / n

print papstats.pformat(t1_2, unit='s', label='t1_2')
print papstats.pformat(t2_2, unit='s', label='t2_2')
print papstats.pformat(T1_2, unit='s', label='T1_2')
print papstats.pformat(T2_2, unit='s', label='T2_2')

print ''
print "## Messung 2"

l = unc.ufloat(0.292, 0.001)
print u"Länge", papstats.pformat(l)

print "Schwingungsperiodendauer"
dt_Schwingung = 0.05
t1_1 = unc.ufloat(39.11, dt_Schwingung) - unc.ufloat(2.20, dt_Schwingung)
t2_1 = unc.ufloat(50.27, dt_Schwingung) - unc.ufloat(13.34, dt_Schwingung)

n = 25

T1_1, T2_1 = t1_1 / n, t2_1 / n

print papstats.pformat(t1_1, unit='s', label='t1_1')
print papstats.pformat(t2_1, unit='s', label='t2_1')
print papstats.pformat(T1_1, unit='s', label='T1_1')
print papstats.pformat(T2_1, unit='s', label='T2_1')

print "Schwebungsperiodendauer"
dt_Schwebung = 0.5
t1_2 = unc.ufloat(64.92, dt_Schwebung) - unc.ufloat(14.06, dt_Schwebung)
t2_2 = unc.ufloat(56.25, dt_Schwebung) - unc.ufloat(5.39, dt_Schwebung)

n = 1.5

T1_2, T2_2 = t1_2 / n, t2_2 / n

print papstats.pformat(t1_2, unit='s', label='t1_2')
print papstats.pformat(t2_2, unit='s', label='t2_2')
print papstats.pformat(T1_2, unit='s', label='T1_2')
print papstats.pformat(T2_2, unit='s', label='T2_2')

print ''
print "## Messung 3"

l = unc.ufloat(0.192, 0.001)
print u"Länge", papstats.pformat(l)

print "Schwingungsperiodendauer"
dt_Schwingung = 0.05
t1_1 = unc.ufloat(26.00, dt_Schwingung) - unc.ufloat(2.29, dt_Schwingung)
t2_1 = unc.ufloat(45.34, dt_Schwingung) - unc.ufloat(21.61, dt_Schwingung)

n = 15

T1_1, T2_1 = t1_1 / n, t2_1 / n

print papstats.pformat(t1_1, unit='s', label='t1_1')
print papstats.pformat(t2_1, unit='s', label='t2_1')
print papstats.pformat(T1_1, unit='s', label='T1_1')
print papstats.pformat(T2_1, unit='s', label='T2_1')

print "Schwebungsperiodendauer"
dt_Schwebung = 0.5
t1_2 = unc.ufloat(76.52, dt_Schwebung) - unc.ufloat(36.71, dt_Schwebung)
t2_2 = unc.ufloat(56.63, dt_Schwebung) - unc.ufloat(17.00, dt_Schwebung)

n = 0.5

T1_2, T2_2 = t1_2 / n, t2_2 / n

print papstats.pformat(t1_2, unit='s', label='t1_2')
print papstats.pformat(t2_2, unit='s', label='t2_2')
print papstats.pformat(T1_2, unit='s', label='T1_2')
print papstats.pformat(T2_2, unit='s', label='T2_2')