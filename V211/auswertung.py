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

def macheSchoenenSchwingungsGraphen(file, xmin, xmax):
    data = np.loadtxt(file + '.txt', skiprows=1)

    title = open(file + '.txt').readline()

    plt.clf()

    x  = data[:, 0]
    y1 = data[:, 1]
    y2 = data[:, 2]

    xrange = (x >= xmin) & (x <= xmax)

    fig = plt.figure()
    ax = fig.add_subplot(111)    # The big subplot
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    ax1.plot(x[xrange], y1[xrange])
    ax2.plot(x[xrange], y2[xrange])

    ax1.yaxis.set_label_position("right")
    ax2.yaxis.set_label_position("right")
    ax1.set_ylabel(u'Pendel 1')
    ax2.set_ylabel(u'Pendel 2')

    ax.set_title(title)
    ax.set_xlabel(ur'Zeit $t$ in $s$')
    ax.set_ylabel(ur'Auslenkung in $^\circ \, (Grad)$')

    ax1.set_xlim([xmin, xmax])

    papstats.savefig_a4(file + '.png')


macheSchoenenSchwingungsGraphen('3.1.1', 0, 10)
macheSchoenenSchwingungsGraphen('3.1.2', 0, 10)
macheSchoenenSchwingungsGraphen('3.2.1', 0, 10)
macheSchoenenSchwingungsGraphen('3.2.2', 0, 10)
macheSchoenenSchwingungsGraphen('3.3.1', 0, 10)
macheSchoenenSchwingungsGraphen('3.3.2', 0, 10)

macheSchoenenSchwingungsGraphen('4.1', 0, 50)
macheSchoenenSchwingungsGraphen('4.2', 0, 30)
macheSchoenenSchwingungsGraphen('4.3', 0, 15)