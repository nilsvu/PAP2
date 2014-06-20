# -*- coding: utf-8 -*-

'''
Auswertung des Versuchs 233 - Fourieroptik
Physikalisches Anfängerpraktikum II, Universität Heidelberg

Autor: Nils Fischer
'''

import numpy as np
import scipy.constants as const
import scipy.integrate as itg
import matplotlib.pyplot as plt
import uncertainties as unc
import uncertainties.unumpy as unp

import os
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0,parentdir)
import papstats


# Eichung

d = unp.uarray([2. * 71.3, 2. * 52.8, 2. * 31.0], 2 * 0.2) * const.milli / 100.
d_px = unp.uarray([1380.44, 1297.31, 1188.49], 5) - unp.uarray([710.88, 803.07, 908.87], 5)
px = np.mean(d / d_px)
print "Eichung: " + papstats.pformat(px / const.milli, unit='mm/px') + " <=> " + papstats.pformat(1. / px * const.milli, unit='px/mm')
px_erw = 14 * const.micro

#####
print "# 3.1: Quantitative Beobachtungen am Einzelspalt"
#####

# Messung bei hoher Intensität
n, x, dx, I, dI = np.loadtxt('2.1.2.txt', skiprows=1, unpack=True)
x = unp.uarray(x, dx)
I = unp.uarray(I, dI)
# Messung bei niedriger Intensität
n_low, x_low, dx_low, I_low, dI_low = np.loadtxt('2.1.1.txt', skiprows=1, unpack=True)
x_low = unp.uarray(x_low, dx_low)
I_low = unp.uarray(I_low, dI_low)

def i_center(x):
	return (len(x) - 1) / 2

# Abstände vom 0. Maximum
i_0 = i_center(n)
i_0_low = i_center(n_low)
x_0 = x_low[i_0_low]
x = x - x_0
x[i_0] = unc.ufloat(0, 0)
x_low = x_low - x_0

# Untergrund
I = I - np.mean(unp.uarray([223.52, 220.72], 5))
I_low = I_low - np.mean(unp.uarray([156.4, 128.21], 5))

def normiere(I):
	return I / I[i_center(I)]

# Normierung
I_low = normiere(I_low)
I = I / np.mean([I[i_0+2], I[i_0-2]]) * np.mean([I_low[i_0_low+2], I_low[i_0_low-2]]) / I_low[i_0_low]
I[i_0] = I_low[i_0_low]

# Merge into combined data
for i_low in range(len(I_low)):
	i = i_0 - i_0_low + i_low
	x[i] = np.mean([x[i], x_low[i_low]])
	I[i] = np.mean([I[i], I_low[i_low]])

# Ordnung
def nspace_centered(x):
	i_0 = i_center(x)
	n = (np.arange(len(x)) - i_0) / 2.
	return n

def nspace_centered_offset(x):
	n = nspace_centered(x)
	for i in range(len(n)):
		if i < i_0:
			n[i] = n[i] - 0.5
		elif i > i_0:
			n[i] = n[i] + 0.5
	return n
n = nspace_centered_offset(I)

print papstats.table(labels=['n', 'x', 'I'], units=[None, 'px', None], columns=[n, x, I])

def fit_linear(x, m):
	return m * x

def approximate_I(n):
	x = n*const.pi + 1e-6
	return np.sin(x)**2/x**2


# Fit Minima
popt, pstats = papstats.curve_fit(fit_linear, n[1::2], x[1::2])

# Plot Abstand
plt.clf()
plt.title(u'Diagramm 3.1: Abstand der Interferenzmaxima und -minima vom Hauptmaximum')
papstats.plot_data(n[1::2], x[1::2], c='b', label='Minima')
papstats.plot_data(n[0::2], x[0::2], c='r', label='Maxima')
papstats.plot_fit(fit_linear, popt, xspace=n, eq='x_{Min}=m*n', punits=['px'])
plt.legend(loc='upper left')
plt.xlabel('Ordnung $n$')
plt.ylabel('Abstand vom Hauptmaximum $x \, [px]$')
papstats.savefig_a4('3.1.png')

# Spaltweite
d = 80 * const.milli * 635 * const.nano / (popt[0] * px) # TODO: use px_erw?
print papstats.pformat(d / const.micro, label='Spaltweite d', unit='um')

# Maxima nach Fit
n_maxima_fit = x[0::2] / popt[0]
diff = np.abs(n[0::2]-n_maxima_fit)
print papstats.table(labels=['n_erw', 'n', u'|n-n_erw|', u'σ-Bereich'], columns=[n[0::2], n_maxima_fit, diff, unp.nominal_values(diff)/unp.std_devs(diff)])

# Plot Intensität
plt.clf()
plt.title(u'Diagramm 3.2: Normierte Intensitätsmaxima und -minima')
papstats.plot_data(n, I, label='Messwerte')
nspace = np.linspace(n[0], n[-1], 300)
Ispace = approximate_I(nspace)
plt.plot(nspace, Ispace, c='black', label='Erwartungswerte')
plt.legend()
plt.xlabel('Ordnung $n$')
plt.ylabel(ur'Normierte Intensität $\frac{I}{I_0}$')
papstats.savefig_a4('3.2.png')

# Vergleich der Intensitäten mit Erwartungswert
def compare_I_max(I, I_erw, nspace):
	I = normiere(I)
	I_max = I[0::2]
	n = nspace(I)
	I_max_erw = I_erw(n[0::2])
	diff = np.abs(I_max-I_max_erw)
	print papstats.table(labels=['Ordnung', 'I', 'I_erw', u'|I-I_erw|', u'σ-Bereich'], columns=[np.arange(len(I_max))-(len(I_max)-1)/2., I_max, I_max_erw, diff, unp.nominal_values(diff)/unp.std_devs(diff)])
compare_I_max(I, approximate_I, nspace_centered_offset)

#####
print "# 3.2: Beugungsstruktur des Doppelspalts"
#####
n, x, dx, I, dI = np.loadtxt('2.2.txt', skiprows=1, unpack=True)
x = unp.uarray(x, dx)
I = unp.uarray(I, dI) - np.mean(unp.uarray([131.58, 134.58], 5))

# Verhältnis g/d der optischen Abbildung
d_DS = unc.ufloat(160.34, 10 / np.sqrt(2)) * px_erw
g_DS = unc.ufloat(224.47, 15) * px_erw + d_DS
V_DS = g_DS/d_DS
print papstats.pformat(V_DS, label=u'Verhältnis g/d')


def approximate_I_doppelspalt(n):
	global V
	return np.cos(n * const.pi)**2 * approximate_I(n / V)	

# Plot Doppelspaltstruktur
plt.clf()
plt.title('Diagramm 3.3: Erwartete Beugungsstruktur am Doppelspalt')
nspace = np.linspace(-10, 10, 1000)
V = V_DS.nominal_value
Ispace_doppelspalt = approximate_I_doppelspalt(nspace)
Ispace_spalt = approximate_I(nspace / V)
V = V_DS.nominal_value + V_DS.std_dev
Ispace_doppelspalt_high = approximate_I_doppelspalt(nspace)
V = V_DS.nominal_value - V_DS.std_dev
Ispace_doppelspalt_low = approximate_I_doppelspalt(nspace)
plt.plot(nspace, Ispace_doppelspalt, c='black', label='Doppelspaltstruktur')
plt.fill_between(nspace, Ispace_doppelspalt_high, Ispace_doppelspalt_low, color='blue', alpha=0.3, lw=0)
plt.plot(nspace, Ispace_spalt, c='red', label=u'Einhüllende Einzelspaltstruktur')
nc = nspace_centered(I)
papstats.plot_data(nc, normiere(I), c='blue', label='Messwerte')
plt.xlabel('Ordnung der Doppelspaltstruktur $n$')
plt.ylabel(ur'relative Intensität $\frac{I}{I_0}$')
plt.legend()
papstats.savefig_a4('3.3.png')

# Vergleiche Intensitäten
	
compare_I_max(I, approximate_I_doppelspalt, nspace_centered)


#####
print "# 3.3: Objektbild als Fouriersynthese"
#####

dx_max = np.loadtxt('2.3.1.txt', skiprows=1, unpack=True)
dx_min = np.loadtxt('2.3.2.txt', skiprows=1, unpack=True)

dx = np.hstack([dx_max, dx_min])

# Äquidistanz
dx_mean = unc.ufloat(np.mean(dx), np.std(dx))
print papstats.pformat(dx_mean, label='MW dx', unit='px')
print papstats.table(labels=['Abstand', 'Abweichung vom MW', u'σ-Bereich'], columns=[dx, dx - dx_mean, unp.nominal_values(np.abs(dx - dx_mean)) / dx_mean.std_dev])

# Spaltweite
b = unc.ufloat(80.7, 0.5) * const.centi
f = unc.ufloat(80, 2) * const.milli
B = b / f - 1
print papstats.pformat(B, label=u'Abbildungsmaßstab B')
d_FS = unc.ufloat(160.72, 5) * px_erw / B
print papstats.pformat(d_FS / const.micro, label=u'Spaltweite d', unit='um')
papstats.print_rdiff(d_FS, d)

# 
def F_einzelspalt(n):
	global y
	return np.sin(const.pi * n) / (const.pi * n) * np.cos(2 * const.pi * n * y)

plt.clf()
plt.suptitle(u'Diagramm 3.4: Objektbild als Fouriersynthese bis zur $n$-ten Ordnung')
yspace = np.linspace(-1, 1, 300)
ax = None
I_0 = None
for n in range(3):
	ax_n = plt.subplot(1, 3, n + 1, sharey=ax)
	plt.title('$n='+str(n+1)+'$')
	if ax is None:
		ax = ax_n
	Ispace = []
	for y in yspace:
		a = itg.quad(F_einzelspalt, 0, n + 1)
		Ispace.append(a[0]**2)
	if I_0 is None:
		I_0 = np.max(Ispace)
	Ispace = Ispace / I_0
	plt.plot(yspace, Ispace)
	plt.xlabel(r'$\frac{y}{d}$')
	if n == 0:
		plt.ylabel(r'$\frac{I}{I_0}$')
papstats.savefig_a4('3.4.png')

# Vergleiche Intensitäten
I_0 = unc.ufloat(3865.11, 35) - np.mean(unp.uarray([302.18, 304.15], 5))
I_1 = unp.uarray([3556.07, 2508.52, 3526.58], [50, 25, 30]) - np.mean(unp.uarray([136.76, 131.93], 5))
I_2 = unp.uarray([3540.55, 2626.71, 3248.80, 2638.55, 3440.44], [30, 20, 50, 10, 65]) - np.mean(unp.uarray([126.30, 133.27], 5))

I = np.hstack([I_0, I_1, I_2]) / I_0
#I_erw = unp.uarray([1, 1, 0.65, 1, 1, 0.7, 0.92, 0.7, 1], 0.01)
I_erw = unp.uarray([1, 0.9, 0.59, 0.9, 0.89, 0.61, 0.81, 0.61, 0.89], 0.02)

diff = np.abs(I - I_erw)
print papstats.table(labels=['I', 'I_erw', u'|I-I_erw|', u'σ-Bereich'], columns=[I, I_erw, diff, unp.nominal_values(diff)/unp.std_devs(diff)])


#####
print "# 3.4: Fourierbild des Doppelspalts"
#####

print u"Fall a: Zulassen der ersten Beugungsmaxima"

d_DS = d_DS / B
g_DS = g_DS / B
print papstats.pformat(d_DS / const.micro, label='Spaltweite d_DS', unit='um')
print papstats.pformat(g_DS / const.micro, label='Spaltabstant g_DS', unit='um')
k_y_1 = 2*const.pi / d_DS
print papstats.pformat(k_y_1, label='Integrationsgrenze k_y_1')

def F_doppelspalt(k_y):
	global y, d, g
	return 2 * d / const.pi * np.cos(k_y * g / 2) * np.sin(k_y * d / 2) / (k_y * d / 2) * np.cos(k_y * y)

yspace = np.linspace(-g_DS.nominal_value, g_DS.nominal_value, 300)
d = d_DS.nominal_value
g = g_DS.nominal_value
Ispace = []
for y in yspace:
	a = itg.quad(F_doppelspalt, 0, k_y_1.nominal_value)
	Ispace.append(a[0]**2)

plt.clf()
plt.title('Diagramm 3.5: Objektbild des Doppelspalts bei Fall a)')
plt.plot(yspace / g * 2, Ispace / np.max(Ispace))
plt.xlabel(r'$y*\frac{2}{g}$')
plt.ylabel(r'$\frac{I}{I_0}$')
papstats.savefig_a4('3.5.png')

k_y_1_erw = 2 * const.pi / (635 * const.nano) * unc.ufloat(22, 0.2) * const.milli / 100 / f
papstats.print_rdiff(k_y_1_erw, k_y_1)

print u"Fall b: Doppelstruktur verschwindet"


k_y_2_erw = 2 * const.pi / (635 * const.nano) * unc.ufloat(6.7, 0.2) * const.milli / 100 / f
k_y_2_begin = k_y_1.nominal_value
k_y_2_end = 1e3
i = 0
plt.clf()
I_0 = None
ax = None
kspace = np.linspace(k_y_2_begin, k_y_2_end, 9)
for k_y in kspace:
	Ispace = []
	for y in yspace:
		a = itg.quad(F_doppelspalt, 0, k_y)
		Ispace.append(a[0]**2)
	ax_i = plt.subplot(3, 3, i + 1, sharey=ax, sharex=ax)
	plt.title('$'+papstats.pformat(k_y*const.milli, label='k_y', unit='mm^{-1}', format='l')+'$')
	if ax is None:
		ax = ax_i
	if I_0 is None:
		I_0 = np.max(Ispace)
	plt.plot(yspace / g * 2, Ispace / I_0)
	if i < 6:
		plt.setp(ax_i.get_xticklabels(), visible=False)
	else:
		plt.xlabel(r'$y*\frac{2}{g}$')
	if i % 3 != 0:
		plt.setp(ax_i.get_yticklabels(), visible=False)
	else:
		plt.ylabel(r'$ \frac{I}{I_0}$')
	i = i+1
papstats.savefig_a4('3.6.png')

k_y_2 = unc.ufloat(kspace[6], (k_y_2_begin-k_y_2_end)/9)
print papstats.pformat(k_y_2, label='Integrationsgrenze k_y_2')
papstats.print_rdiff(k_y_2_erw, k_y_2)




