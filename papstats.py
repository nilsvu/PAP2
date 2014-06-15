# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats as st
import scipy.optimize as opt
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import prettytable as pt
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

def plot_data(xdata, ydata, ax=plt, **kwargs):
    xerr = unp.std_devs(xdata)
    if np.sum(xerr)==0:
        xerr = None
    yerr = unp.std_devs(ydata)
    if np.sum(yerr)==0:
        yerr = None
    if not (kwargs.has_key('ls') or kwargs.has_key('linestyle')):
        kwargs['ls'] = 'none'
    if not kwargs.has_key('marker'):
        kwargs['marker'] = '.'
    return ax.errorbar(unp.nominal_values(xdata), unp.nominal_values(ydata), xerr=xerr, yerr=yerr, **kwargs)


def plot_fit(fit, popt, pstats=None, xspace=np.linspace(0, 1), xscale=1., yscale=1., eq=None, plabels=None, punits=None, **kwargs):
    if eq is None:
        eq = ''
    else:
        if eq.find('$') == -1:
            eq = '$' + eq + '$ '
    if plabels is None:
        plabels = inspect.getargspec(fit)[0][1:]
    if punits is None:
        punits = [None for plabel in plabels]
    if not kwargs.has_key('label'):
        kwargs['label'] = 'Fit ' + eq + 'mit:' + ''.join(['\n$%s$' % pformat(popt[i], label=plabels[i], unit=punits[i], format='l') for i in range(len(popt))]) + ('\n' + pstats.legendstring() if pstats is not None else '')
    ydata = fit(xspace, *unp.nominal_values(popt))
    return plt.plot(xspace * xscale, ydata * yscale, **kwargs)


def savefig_a4(filename):
    fig = plt.gcf()
    fig.set_size_inches(11.69, 8.27)
    plt.savefig(filename, dpi=150)


# (Pretty)Tables

def table(labels=None, units=None, columns=None, filename=None, rowlabels=None, **kwargs):
    table = pt.PrettyTable()
    if rowlabels is not None:
        table.add_column('', rowlabels)
    for i in range(len(columns)):
        try:
            label = labels[i]
        except TypeError:
            label = None
        try:
            unit = units[i]
        except TypeError:
            unit = None
        if label is not None and unit is not None:
            label += ' ['+unit+']'
        table.add_column(label, pformat(columns[i], format='p', **kwargs))
    if filename is not None:
        with open(filename, 'w') as file:
            file.write(table.get_string())
    return table


# Formatting

def pformat(v, dv=None, prec=2, label=None, unit=None, format=None, tex=False):
    if isinstance(v, basestring):
        return v
    try:
        return np.array([pformat(vi, dv=dv, prec=prec, label=label, unit=unit, format=format) for vi in v])
    except TypeError:
        pass

    if label is None and isinstance(v, unc.Variable):
        label = v.tag
    if label is None:
        label = ''
    else:
        label += '='
    if unit is None:
        unit = ''

    # use uncertainties module formatting
    if isinstance(v, unc.UFloat):
        if format is None or format=='c' or format=='console' or format=='p' or format=='pretty':
            format = '.' + str(prec) + 'uP'
        if format=='l' or format=='latex':
            format = '.' + str(prec) + 'uL'
            
    if format is None or format=='l' or format=='latex' or format=='c' or format=='console' or format=='p' or format=='pretty':
        format = '.' + str(prec) + 'e'

    return label + (u'{:' + format + u'}').format(v) + unit

    # format numbers without uncertainties

    if unit is None:
        unit = ''
    if label is None:
        label = ''
    else:
        label += '='

    e = np.floor(np.log10(np.abs(v)))
    o = 10 ** (e - prec + 1)
    v = round_ordnung(v, o)
    if dv is not None:
        dv = round_ordnung(dv, o)
        return (ur"%s(%." + str(prec - 1) + "f\pm%." + str(prec - 1) + "f)*10^{%d}%s") % (label, v / 10 ** e, dv / 10 ** e, e, unit)
    else:
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

def rdiff(r, r_erw=0):
    d = np.abs(r - r_erw)
#    if not isinstance(r_erw, unc.UFloat):
#        r_erw = unc.ufloat(r_erw, 0)
#    return d, (d / (r_erw if r_erw.n!=0 else r) if np.abs(r_erw.n)+np.abs(r.n)!=0 else None) , d.n / d.s
    return d, d / r_erw, unp.nominal_values(d) / unp.std_devs(d)

def print_rdiff(r, r_erw):
    d, drel, s = rdiff(r, r_erw)
    print 'Berechnung: ' + pformat(r)
    print 'Erwartungswert: ' + pformat(r_erw)
    print 'Abweichung: ' + pformat(d)
    print 'rel. Abweichung: ' + pformat(drel * 100, unit='%')
    print 'Sigmabereich: ' + pformat(s)