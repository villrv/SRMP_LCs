import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize

time,flux = np.loadtxt('fake_data.dat',usecols=(0,1),unpack=True)

def func(x, m, b):
		return m*x + b

popt, pcov = curve_fit(func, time, flux)

sigs = np.zeros(len(time))+0.1


def chi2(perams, time, flux):
		m = perams[0]
		b = perams[1]
		x = time
		y = flux
		return sum(((y-func(x, m, b))**2)/(sigs**2))

x0 = [.5, 2]

res = minimize(chi2, x0, args=(time,flux))

plt.plot(time, func(time, res.x[0], res.x[1]))
plt.plot(time,flux,'o')
plt.show()

