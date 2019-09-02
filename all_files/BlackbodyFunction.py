import numpy as np


def blackbody(r, temp, wavelength):
	sB = 5.670*(10*8)
	c = 3.0*(10**10)
	h = 6.626*(10**-27)
	k = 1.38*(10**-16)

	planks_function = ((2*np.pi*(c**2)*h)/(wavelength**5))*(1/((np.exp((h*c)/(wavelength*k*temp)))-1))

	y = ((2*np.pi*(c**2)*h)/(wavelength**5))*(1/((np.exp((h*c)/(wavelength*k*temp)))-1))

	flux = y*(r**2)
	return flux 