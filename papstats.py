import numpy as np
import scipy.stats as st
import scipy.optimize as opt
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import re

# Fit

def curve_fit(fit, xdata, ydata, p0=None):
    popt, pcov = opt.curve_fit(fit, unp.nominal_values(xdata), unp.nominal_values(ydata), sigma=unp.std_devs(xdata+ydata), p0=p0)
    popt = unp.uarray(popt, np.sqrt(np.diagonal(pcov)))
    pstats = PAPStats(unp.nominal_values(ydata), fit(unp.nominal_values(xdata), *unp.nominal_values(popt)), sigma=unp.std_devs(xdata+ydata), ddof=len(popt))
    return popt, pstats

def chisquared(ydata, ymodel, sigma=1, ddof=0):
    chisq = np.sum(((ydata-ymodel)/sigma)**2)
    return chisq, chisq/(ydata.size-ddof)

# Statistics

class PAPStats:
    def __init__(self, ydata, ymodel, sigma=1, ddof=0):
        self.chisq = chisquared(ydata, ymodel, sigma, ddof)
        self.pearsonr = st.pearsonr(ydata, ymodel)
        self.rsquared = self.pearsonr[0]**2
        self.residue = ymodel-ydata
    
    def __str__(self):
        return "< Chi-Squared: "+str(self.chisq)+", Pearson R: "+str(self.pearsonr)+", R Squared: "+str(self.rsquared)+" >"
    
    def __repr__(self):
        return str(self)
    
    def legendstring(self):
        return '$\chi^2=%.2f$, $\chi^2_{red}=%.2f$, $r^2=%.5f$' % (self.chisq[0], self.chisq[1], self.rsquared)

# Plots

def plot_data(xdata, ydata, **kwargs):
    plt.errorbar(unp.nominal_values(xdata), unp.nominal_values(ydata), xerr=unp.std_devs(xdata), yerr=unp.std_devs(ydata), ls='none', marker='.', **kwargs)

def plot_fit(fit, popt, pstats, xfit, eq=None, **kwargs):
    if eq is None:
        eq = ''
    else:
        eq = '$'+eq+'$ '
    label = 'Fit '+eq+'mit:'+''.join(['\n$%s$' % pformat(p) for p in popt])+'\n'+pstats.legendstring()
    plt.plot(xfit, fit(xfit, *unp.nominal_values(popt)), label=label, **kwargs)

def plot(xdata, ydata, fit=None, popt=None, pstats=None, xfit=None, **kwargs):
    plot_data(xdata, ydata)
    if fit is not None:
        if xfit is None:
            xfit = np.linspace(xdata[0], xdata[-1], num=100)
        plot_fit(fit, popt, pstats, xfit)

def savefig_a4(filename):
    fig = plt.gcf()
    fig.set_size_inches(11.69,8.27)
    plt.savefig(filename, dpi=144)

# Formatting

def pformat(v, dv=None, prec=2, label=None, unit=None):
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
        return label+'{:.2uL}'.format(v)+unit
    
    # format numbers without uncertainties

    if unit is None:
        unit = ''
    if label is None:
        label = ''
    else:
        label += '='

    if dv is not None:
        e = np.floor(np.log10(zahl))
        o = 10**(e-prec+1)
        v = round_ordnung(v, o)
        dv = round_ordnung(dv, o)
        return (r"$%s(%."+str(prec-1)+"f\pm%."+str(prec-1)+"f)*10^{%d}%s$") % (label, v / 10**e, dv / 10**e, e, unit)
    else:
        e = get_e(v)
        o = 10**(e-prec+1)
        v = round_ordnung(v, o)
        string = r"$%s%."+str(prec-1)+"f*10^{%d}%s$"
        return string % (label, v/10.0**e, e,  unit)

def round_ordnung(v, o):
    return np.round(v / o) * o

# Calculations

def print_rdiff(r, r_erw):
    d = np.abs(r-r_erw)
    print 'Berechnung:', r
    print 'Erwartungswert:', r_erw
    print 'Abweichung: {:.2u}'.format(d)
    print 'rel. Abweichung: {:.2}'.format(d.n/r_erw.n*100)
    print 'Sigmabereich: {:.2}'.format(d.n/d.s)
