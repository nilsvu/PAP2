# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 256 - Röntgenfluoreszenz
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import scipy.constants as const
import uncertainties as unc
import uncertainties.unumpy as unp

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats


#####
print u"\n# 3.1: Betrieb als Kältemaschine"
#####

# Kompensations-Heizleistung entspricht der Kälteleistung
U_H = unc.ufloat(4.32, 0.01)
I_H = 5 * unc.ufloat(0.91, 0.01)
P_K = U_H * I_H
print papstats.pformat(P_K, format='c', unit='W', label=u'Kälteleistung der Kältemaschine P_K')

print u"Energiebilanz:"

# zugeführte Kompensations-Wärmemenge pro Umdrehung
f_M = unc.ufloat(287.6, 1) / const.minute
Q_2 = P_K / f_M
print papstats.pformat(Q_2, format='c', unit='J', label='Q_2')

# an Kühlwasser abgegebene Wärmemenge
c_W = 4180
rho_W = 998.2
dT = unc.ufloat(2.5, 0.1)
V = unc.ufloat(240, 2) * const.milli * const.liter / const.minute
Q_1 = c_W * rho_W * dT * V / f_M
print papstats.pformat(Q_1, format='c', unit='J', label='Q_1')

print papstats.pformat(Q_1 - Q_2, format='c', unit='J', label=u'zugeführte Energie nach Differenz W')

# zugeführte Motorleistung der Kältemaschine pro Umdrehung nach Messung 2.2
U_M = unc.ufloat(24.0, 0.1)
I_M = unc.ufloat(1.5, 0.1)
f_M = unc.ufloat(281.6, 0.1) / const.minute
W_M = U_M * I_M / f_M
print papstats.pformat(W_M, format='c', unit='J', label=u'tatsächliche Motorleistung W_M')


#####
print u"\n# 3.2: Betrieb als Kältemaschine und Wärmepumpe"
#####

## 3.2.1: Kältemaschine

# latente Wärme
m_W = unc.ufloat(1, 0.5) * const.gram
Q_l = 335 / const.gram * m_W
# Gefrierzeit
dt = unc.ufloat(400, 10) - unc.ufloat(180, 10)
print papstats.pformat(dt, format='c', unit='s', label='Gefrierzeit dt')
# Kälteleistung der Kältemaschine
P_K = Q_l / dt
print papstats.pformat(P_K, format='c', unit='W', label=u'Kälteleistung der Kältemaschine P_K')


#####
print u"\nBetrieb als Wärmekraftmaschine"
#####

f = np.mean(unp.uarray([340.3, 340.4, 340.7], [0.5, 0.3, 0.4])) / const.minute
print papstats.pformat(f, format='c', unit='s^(-1)', label='f')
P_el = unc.ufloat(13.21, 0.01) * unc.ufloat(2.9, 0.01) * 5
print papstats.pformat(P_el, format='c', unit='W', label='P_el')
Q_el = P_el / f_M
print papstats.pformat(Q_el, format='c', unit='J', label='Q_el')
P_ab = c_W * rho_W * (unc.ufloat(24.1, 0.2) - unc.ufloat(17.7, 0.1)) * (unc.ufloat(240, 2) * const.milli * const.liter / const.minute)
print papstats.pformat(P_ab, format='c', unit='W', label='P_ab')
Q_ab = P_ab / f_M
print papstats.pformat(Q_ab, format='c', unit='J', label='Q_ab')
Q_pV = np.mean(np.array([26503, 26728, 26781])) * 1e-4
print papstats.pformat(Q_pV, format='c', unit='J', label='Q_pV')
P_pV = Q_pV * f_M
print papstats.pformat(P_pV, format='c', unit='W', label='P_pV')
eta = Q_pV / Q_el
print papstats.pformat(eta, format='c', label='eta')

Q_V = Q_el - Q_ab - Q_pV
print papstats.pformat(Q_V, format='c', unit='J', label=u'Motorverluste Q_V')

# TODO: Drehmomentmessung
