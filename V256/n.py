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

Ry = 13.6*const.eV


converters = dict.fromkeys(range(2, 4), uncertainties.ufloat_fromstr)
Z, Ea, Eb = np.loadtxt('elements.txt', skiprows=1, usecols=range(1, 4), unpack=True, converters=converters)

E = np.transpose(np.array(Ea, Eb))*const.kilo*const.eV
ERy = unp.sqrt(E/Ry)

n1 = 1
n2 = np.array([2., 3.])


print np.transpose(Z-np.transpose(ERy/np.sqrt(n1**(-2)-n2**(-2))))




converters = dict.fromkeys(range(3), uncertainties.ufloat_fromstr)
np.loadtxt('legierungen.txt', converters=converters)
