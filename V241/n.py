# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 241 - Wechselstromeigenschaften von RLC-Gliedern
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

