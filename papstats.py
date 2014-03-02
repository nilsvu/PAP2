# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats as st
import scipy.optimize as opt
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import inspect

# Fit

def curve_fit(fit, xdata, ydata, sigma=None, p0=None):
    if sigma is None:
        # TODO: Nur y-Fehler? chisquared vergleicht nur y-Differenzen, x-Fehler relevant für Fit?
        sigma = unp.std_devs(ydata)
        if np.sum(sigma) == 0:
            sigma = None
    xdata = unp.nominal_values(xdata)
    ydata = unp.nominal_values(ydata)
    popt, pcov = opt.curve_fit(fit, xdata, ydata, sigma=sigma, p0=p0)
    popt = unp.uarray(popt, np.sqrt(np.diagonal(pcov)))
    pstats = PAPStats(ydata, fit(xdata, *unp.nominal_values(popt)), sigma=sigma, ddof=len(popt))
    return popt, pstats


def chisquared(ydata, ymodel, sigma=None, ddof=0):
    if sigma is None:
        sigma = 1
    chisq = np.sum(((ydata - ymodel) / sigma) ** 2)
    dof = len(ydata) - ddof
    return chisq, chisq / dof, 1 - st.chi2.cdf(chisq, dof)


# Statistics

class PAPStats:
    def __init__(self, ydata, ymodel, sigma=None, ddof=0):
        self.chisq = chisquared(ydata, ymodel, sigma, ddof)
        self.pearsonr = st.pearsonr(ydata, ymodel)
        self.rsquared = self.pearsonr[0] ** 2
        self.residual = ymodel - ydata

    def __str__(self):
        return "< Chi-Squared: " + str(self.chisq) + ", Pearson R: " + str(self.pearsonr) + ", R Squared: " + str(
            self.rsquared) + " >"

    def __repr__(self):
        return str(self)

    def legendstring(self):
        return '$\chi^2=%.2f$, $\chi^2_{red}=%.2f$, $p=%.2f' % (
        self.chisq[0], self.chisq[1], self.chisq[2] * 100) + r'\%$'


# Plots

def plot_data(xdata, ydata, **kwargs):
    plt.errorbar(unp.nominal_values(xdata), unp.nominal_values(ydata), xerr=unp.std_devs(xdata),
                 yerr=unp.std_devs(ydata), ls='none', marker='.', **kwargs)


def plot_fit(fit, popt, pstats, xspace, xscale=1., yscale=1., eq=None, plabels=None, punits=None, **kwargs):
    if eq is None:
        eq = ''
    else:
        eq = '$' + eq + '$ '
    if plabels is None:
        plabels = inspect.getargspec(fit)[0][1:]
    if punits is None:
        punits = [None for plabel in plabels]
    label = 'Fit ' + eq + 'mit:' + ''.join(['\n$%s$' % pformat(popt[i], label=plabels[i], unit=punits[i]) for i in range(len(popt))]) + '\n' + pstats.legendstring()
    ydata = fit(xspace, *unp.nominal_values(popt))
    plt.plot(xspace * xscale, ydata * yscale, label=label, **kwargs)


def savefig_a4(filename):
    fig = plt.gcf()
    fig.set_size_inches(11.69, 8.27)
    plt.savefig(filename, dpi=150)


# Formatting

def pformat(v, dv=None, prec=2, label=None, unit=None, format=None):
    # use uncertainties module formatting
    if isinstance(v, unc.UFloat):
        if label is None and isinstance(v, unc.Variable):
            label = v.tag
        if label is None:
            label = ''
        else:
            label += '='
        if unit is None:
            unit = ''
        if format is None:
            format = '.' + str(prec) + 'uL'
        return label + ('{:' + format + '}').format(v) + unit

    # format numbers without uncertainties

    if unit is None:
        unit = ''
    if label is None:
        label = ''
    else:
        label += '='

    if dv is not None:
        e = np.floor(np.log10(v))
        o = 10 ** (e - prec + 1)
        v = round_ordnung(v, o)
        dv = round_ordnung(dv, o)
        return (ur"%s(%." + str(prec - 1) + "f\pm%." + str(prec - 1) + "f)*10^{%d}%s") % (
        label, v / 10 ** e, dv / 10 ** e, e, unit)
    else:
        e = np.floor(np.log10(v))
        o = 10 ** (e - prec + 1)
        v = round_ordnung(v, o)
        if np.abs(e) > prec:
            string = r"%s%." + str(prec - 1) + r"f \times 10^{%d}%s"
            string = string % (label, v / 10.0 ** e, e, unit)
        else:
            string = r"%s%." + str(prec) + "g%s"
            string = string % (label, v, unit)
        return string


def round_ordnung(v, o):
    return np.round(v / o) * o


# Calculations

def rdiff(r, r_erw):
    d = np.abs(r - r_erw)
    return d, d / r_erw, d.n / d.s

def print_rdiff(r, r_erw):
    d, drel, s = rdiff(r, r_erw)
    print 'Berechnung: {:.2u}'.format(r)
    print 'Erwartungswert: {:.2u}'.format(r_erw)
    print 'Abweichung: {:.2u}'.format(d)
    print 'rel. Abweichung: {:.2u}%'.format(drel * 100)
    print 'Sigmabereich: {:.2}'.format(s)