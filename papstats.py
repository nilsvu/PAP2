import numpy as np
import scipy.stats as st

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

def pformat(v, dv=None, label=None, signi=2, unit=None):
    if unit is None:
        unit = ''
    if label is None:
        label = ''
    else:
        label += '='
    if dv is not None:
        e = get_e(dv)
        o = 10**(e-signi+1)
        v = round_ordnung(v, o)
        dv = round_ordnung(dv, o)
        return (r"$%s(%."+str(signi-1)+"f\pm%."+str(signi-1)+"f)*10^{%d}%s$") % (label, v / 10**e, dv / 10**e, e, unit)
    else:
        e = get_e(v)
        o = 10**(e-signi+1)
        v = round_ordnung(v, o)
        string = r"$%s%.1f*10^{%d}%s$"
        return string % (label, v/10.0**e, e,  unit)

def round_ordnung(v, o):
    return np.round(v / o) * o

def get_e(zahl):
    return np.floor(np.log10(zahl))