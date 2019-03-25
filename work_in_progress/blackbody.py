import numpy as np
import matplotlib.pyplot as plt

σB = 5.670*(10*8)
c = 299792458
h = 6.62607004*(10*34)
k = 1.38064852*(10*23)
wavelength = 4*(10**-5)
T = 5000


planks_function = ((2*np.pi*(c**2)*h)/(wavelength**5))*(1/((np.exp((h*c)/(wavelength*k*T)))-1))

x = linspace(3*(10**-5), 5*(10**-5), 100)
y = ((2*np.pi*(c**2)*h)/(x**5))*(1/((np.exp((h*c)/(x*k*T)))-1))

plot(x,y)