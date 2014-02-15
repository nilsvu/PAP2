import numpy as np

def chisquared(ydata, ymodel, std=1):
    return np.sum(((ydata-ymodel)/std)**2)

def chisquared_red(ydata, ymodel, ddof=2, std=1):
    return chisquared(ydata, ymodel, std)/(ydata.size-1-ddof)