import numpy as np
import scipy.stats as st
import uncertainties as unc

def chisquared(ydata, ymodel, sigma=1, ddof=0):
    chisq = np.sum(((ydata-ymodel)/sigma)**2)
    return chisq, chisq/(ydata.size-ddof)

class PAPStats:
    def __init__(self, ydata, ymodel, sigma=1, ddof=0):
        self.chisq = chisquared(ydata, ymodel, sigma, ddof)
        self.pearsonr = st.pearsonr(ydata, ymodel)
        self.rsquared = self.pearsonr[0]**2

    def __str__(self):
        return "< Chi-Squared: "+str(self.chisq)+", Pearson R: "+str(self.pearsonr)+", R Squared: "+str(self.rsquared)+" >"

    def __repr__(self):
        return str(self)

    def legendstring(self):
        return '$\chi^2=%.2f$, $\chi^2_{red}=%.2f$, $r^2=%.5f$' % (self.chisq[0], self.chisq[1], self.rsquared)


def pformat(v, dv=None, prec=2, label=None, unit=None):
    if unit is None:
        unit = ''
    if label is None:
        label = ''
    else:
        label += '='

    # use uncertainties module formatting
    if isinstance(v, unc.Variable) or isinstance(v, unc.AffineScalarFunc):
        return '$'+label+ '{:.2uL}'.format(v) +unit+'$'
    
    # format numbers without uncertainties
    if dv is not None:
        e = get_e(dv)
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

def get_e(zahl):
    return np.floor(np.log10(zahl))

def print_rdiff(r, r_lit):
    d = np.abs(r-r_lit)
    print 'Abweichung: {:.2u}\nrel. Abweichung: {:.2}%\nSigmabereich: {:.2}'.format(d, d.n/r_lit.n*100, d.n/d.s)
