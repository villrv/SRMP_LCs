import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize

time,flux = np.loadtxt('fake_data.dat',usecols=(0,1),unpack=True)

sigs = np.zeros(len(time))+0.1

def func(x, m, b):
	return m * x + b

def chi2(params, time, flux):
	m = params[0]
	b = params[1]

	# args[0] = x
	# args[1] = y

	return np.sum(((flux - func(time, m, b)**2)/sigs))

x0 = [2, 0.5]

res = minimize(chi2, x0, args=(time, flux))

# popt, pcov = curve_fit(func, time, flux)

plt.plot(time, func(time, res.x[0], res.x[1]))
plt.plot(time,flux,'o')
plt.show()

