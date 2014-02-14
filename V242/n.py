import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c
import scipy.stats as st

data = np.loadtxt('1.a.txt')

Ue = data[:,0]
dUe = 2*c.milli

dUrel = 0.03
Ua = data[:,1:3]
dUa = dUrel*Ua

dRrel = 0.05
Rg = np.array([48.7, 274, 680])*c.kilo
dRg = dRrel * Rg

slope, intercept, r_value, p_value, std_err = st.linregress(Ue, Ua[:,0])
slope1, intercept1, r_value1, p_value1, std_err1 = st.linregress(Ue[2:-2], Ua[2:-2,1])

plt.clf()
plt.errorbar(Ue, Ua[:,0], xerr=dUe, yerr=dUa[:,0], ls='none')
plt.plot(Ue, slope*Ue+intercept)
plt.errorbar(Ue,Ua[:,1], xerr=dUe, yerr=dUa[:,1], ls='none')
plt.plot(Ue, slope1*Ue+intercept1)
plt.savefig('1.a.png')