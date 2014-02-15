import numpy as np

def chisquared(ydata, ymodel, std=1, ddof=0):
    chisq = np.sum(((ydata-ymodel)/std)**2)
    return chisq, chisq/(ydata.size-ddof)