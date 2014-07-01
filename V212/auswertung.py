# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 212 - Zähigkeit von Flüssigkeiten
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Autoren: Nils Fischer und Christian Kohlstedde
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
