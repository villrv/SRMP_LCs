import numpy as np
import matplotlib.pyplot as plt

σB = 5.670*(10*8)
c = 3.0*(10**10)
h = 6.626*(10**-27)
k = 1.38*(10**-16)
wavelength = 4*(10**-5)
T = 310


planks_function = ((2*np.pi*(c**2)*h)/(wavelength**5))*(1/((np.exp((h*c)/(wavelength*k*T)))-1))

x = np.linspace(3*(10**-5), 150*(10**-5), 100)
y = ((2*np.pi*(c**2)*h)/(x**5))*(1/((np.exp((h*c)/(x*k*T)))-1))

plt.plot(x,y)

plt.show()