import numpy as np
import scipy.stats as st

def chisquared(ydata, ymodel, std=1, ddof=0):
    chisq = np.sum(((ydata-ymodel)/std)**2)
    return chisq, chisq/(ydata.size-ddof)

class PAPStats:
    def __init__(self, ydata, ymodel, std=1, ddof=0):
        self.chisq = chisquared(ydata, ymodel, std, ddof)
        self.pearsonr = st.pearsonr(ydata, ymodel)
        self.rsquared = self.pearsonr[0]**2

    def __str__(self):
        return "< Chi-Squared: "+str(self.chisq)+", Pearson R: "+str(self.pearsonr)+", R Squared: "+str(self.rsquared)+" >"

    def __repr__(self):
        return str(self)

    def legendstring(self):
        return '$\chi^2=%.2f$, $\chi^2_{red}=%.2f$, $r^2=%.5f$' % (self.chisq[0], self.chisq[1], self.rsquared)