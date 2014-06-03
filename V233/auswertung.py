# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 233 - Fourieroptik
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Author: Nils Fischer
'''

import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import uncertainties as unc
import uncertainties.unumpy as unp

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats


# Eichung

d = unp.uarray([2. * 71.3, 2. * 52.8, 2. * 31.0], np.sqrt(2) * 0.2) * const.milli / 100.
d_px = unp.uarray([1380.44, 1297.31, 1188.49], 5) - unp.uarray([710.88, 803.07, 908.87], 5)
px = np.mean(d / d_px)
print "Eichung: " + papstats.pformat(px / const.milli, unit='mm/px') + " <=> " + papstats.pformat(1. / px * const.milli, unit='px/mm')


#####
print "# 3.1: Quantitative Beobachtungen am Einzelspalt"
#####

# Messung bei hoher Intensität
o, x, dx, I, dI = np.loadtxt('2.1.2.txt', skiprows=1, unpack=True)
x = unp.uarray(x, dx)
I = unp.uarray(I, dI)

i_0 = (len(o) - 1) / 2 - 1
dx = np.array([np.abs(x[i_0+i] - x[i_0-i]) for i in range(i_0 + 1)])

# Messung bei niedriger Intensität
o_low, x_low, dx_low, I_low, dI_low = np.loadtxt('2.1.1.txt', skiprows=1, unpack=True)
x_low = unp.uarray(x_low, dx_low)
I_low = unp.uarray(I_low, dI_low)

# Normierungsfaktoren
I_0_low = I_low[4] # 0. Maximum bei niedriger Intensität
I_0 = I_0_low / np.mean([I_low[2], I_low[6]]) * np.mean()

for i_low in range(len(o_low)):
    i_high = i - 3
    if o_low[i_low] == 0:
        x[i_high] = x_low[i_low]
    




