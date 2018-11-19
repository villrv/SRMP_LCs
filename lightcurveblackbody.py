''' Present an interactive function explorer with slider widgets.

Scrub the sliders to change the properties of the ``sin`` curve, or
type into the title text box to update the title of the plot.

Use the ``bokeh serve`` command to run the example by executing:

    bokeh serve sliders.py

at your command prompt. Then navigate to the URL

    http://localhost:5006/sliders

in your browser.

'''
import scipy.integrate as integrate
import numpy as np
import requests

from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput
from bokeh.plotting import figure, output_file, show

output_file("line.html")

# Set up data
def querry_single_osc(object_name):
    '''
    Download a light curve from the Open Supernova Catalog.
    The output will be in MJD and Magnitude
    '''
    osc_link    = 'https://astrocats.space/api/' + object_name + '/photometry/time+magnitude+e_magnitude+band+upperlimit'
    osc_request = requests.get(osc_link).json()
    osc_data    = np.array(osc_request[object_name]['photometry'])
    photometry_time  = osc_data.T[0].astype(float)
    photometry_mag   = osc_data.T[1].astype(float)
    photometry_sigma = osc_data.T[2]
    photometry_band  = osc_data.T[3]
    photometry_limit = osc_data.T[4]

    # Convert empty sigmas to -1.0
    photometry_sigma[np.where(photometry_sigma == '')] = -1.0
    photometry_sigma = photometry_sigma.astype(float)

    # Reformat Upper Limit
    detection = np.where(photometry_limit != 'True')
    return photometry_time, photometry_mag, photometry_sigma, photometry_band, detection

def querry_distance(object_name):
    '''
    Download a light curve from the Open Supernova Catalog.
    The output will be in MJD and Magnitude
    '''
    osc_link    = 'https://astrocats.space/api/' + object_name + '/lumdist+redshift'
    osc_request = requests.get(osc_link).json()
    osc_data    = osc_request[object_name]
    lumdist = np.array(osc_data['lumdist'][0]['value']).astype(float)
    redshift = np.array(osc_data['redshift'][0]['value']).astype(float)

    return lumdist, redshift

photometry_time, photometry_mag, photometry_sigma, photometry_band, detection = querry_single_osc("DES17C1ffz")

lumdist, redshift = querry_distance("DES17C1ffz")
print(lumdist)
print(redshift)

N = 200
t_original = np.linspace(0, 100, N)
t = t_original
M = 5*(1.989*(10**33))
f = 0.5
k = 0.25
v = 12*(10**8)
tn = 8.8
B = 13.8
c = 3*(10**10)
E = 3.9*(10**10)
td = (((2*k*M)/(B*c*v))**(1./2.))/60./60./24.
integrand = (t/td)*(np.exp(t**2/td**2))*(np.exp(-t/tn))
my_int = integrate.cumtrapz(t,integrand, initial = 0)
L = ((2*M*f)/(td*24*60*60)) * (np.exp(-t**2)/td**2) * E * my_int


source = ColumnDataSource(data=dict(x=t, y=L, yB=L, yr=L, yi=L, yV=L, yU=L))


# Set up plot
plot = figure(plot_height=400, plot_width=400, title="super cool parabola",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[np.min(photometry_time) - 20, np.max(photometry_time) + 100], y_range=[5, 20])

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
plot.line('x', 'yB', source=source, line_width=3, line_alpha=0.6, color="pink")
plot.line('x', 'yr', source=source, line_width=3, line_alpha=0.6, color="orange")
plot.line('x', 'yi', source=source, line_width=3, line_alpha=0.6, color="blue")
plot.line('x', 'yV', source=source, line_width=3, line_alpha=0.6, color="turquoise")
plot.line('x', 'yU', source=source, line_width=3, line_alpha=0.6, color="purple")

arrayoftimes = np.array(photometry_time)

# Set up widgets
text = TextInput(title="title", value='my parabola')
M_slider = Slider(title="Ejecta Mass", value=5, start=1.0, end=10.0, step=0.5)
f_slider = Slider(title="Fraction of Radio Active Stuff", value=0.0, start=0.0, end=1.0, step=0.1)
v_slider = Slider(title="Ejecta Velocity", value= 12000, start=5000, end=20000, step=1000)
k_slider = Slider(title="Kapa", value=0.25, start= 0.1, end=0.5, step=0.05)
T_slider = Slider(title="Time", value= arrayoftimes.min() - 10, start= arrayoftimes.min() - 10, end=arrayoftimes.max() + 10, step= 10)

count = 0
for x in photometry_time:
	if photometry_band[count] == "r":
		plot.circle(x, photometry_mag[count], size=5, color="orange", alpha=0.5)
	elif photometry_band[count] == "i":
		plot.circle(x, photometry_mag[count], size=5, color="blue", alpha=0.5)
	elif photometry_band[count] == "g":
		plot.circle(x, photometry_mag[count], size=5, color="turquoise", alpha=0.5)
	elif photometry_band[count] == "z":
		plot.circle(x, photometry_mag[count], size=5, color="purple", alpha=0.5)
	elif photometry_band[count] == "B":
		plot.circle(x, photometry_mag[count], size=5, color="pink", alpha=0.5)
	count += 1

print(T_slider.value)

show(plot)


# Set up callbacks
def update_title(attrname, old, new):
    plot.title.text = text.value

text.on_change('value', update_title)

def blackbody(r, temp, wavelength):
	sB = 5.670*(10*8)
	c = 3.0*(10**10)
	h = 6.626*(10**-27)
	k = 1.38*(10**-16)

	planks_function = ((2*np.pi*(c**2)*h)/(wavelength**5))*(1/((np.exp((h*c)/(wavelength*k*temp)))-1))

	y = ((2*np.pi*(c**2)*h)/(wavelength**5))*(1/((np.exp((h*c)/(wavelength*k*temp)))-1))

	flux = y*(r**2)
	return flux 


def update_data(attrname, old, new):

    # Get the current slider values
    M = M_slider.value * 2.e33
    f = f_slider.value
    v = v_slider.value * 1.e5
    k = k_slider.value
    t = t_original
    tn = 8.8
    B = 13.8
    c = 3*(10**10)
    E = 3.9*(10**10)
    td = (((2.*k*M)/(B*c*v))**(1./2.))/60./60./24.
    integrand = (t/td)*(np.exp(t**2/td**2))*(np.exp(-t/tn))
    my_int = integrate.cumtrapz(integrand,t, initial = 0)
    L = ((2.*M*f)/(td)) * (np.exp((-t**2)/td**2)) * E * my_int
    magnitude = -2.5*np.log10(L/4e33)+4.3
    radii = v * t * 24*60*60
    temperature = (L/(4*np.pi*(radii**2)*(5.67*10**-5)))**0.25
    county = 0
    wavelength = 5*10**-5
    wavelength_B = 4.45*10**-5
    wavelength_r = 6.58*10**-5
    wavelength_i = 8.06*10**-5
    wavelength_V = 5.51*10**-5
    wavelength_U = 3.65*10**-5
    distance = lumdist * (3*10**24)
    luminosityblackbody = blackbody(radii, temperature, 5.*(10**-5)*np.ones(len(radii))) * wavelength**2 / c / distance**2
    luminosityblackbody_B = blackbody(radii, temperature, 5.*(10**-5)*np.ones(len(radii))) * wavelength_B**2 / c / distance**2
    luminosityblackbody_r = blackbody(radii, temperature, 5.*(10**-5)*np.ones(len(radii))) * wavelength_r**2 / c / distance**2
    luminosityblackbody_i = blackbody(radii, temperature, 5.*(10**-5)*np.ones(len(radii))) * wavelength_i**2 / c / distance**2
    luminosityblackbody_V= blackbody(radii, temperature, 5.*(10**-5)*np.ones(len(radii))) * wavelength_V**2 / c / distance**2
    luminosityblackbody_U = blackbody(radii, temperature, 5.*(10**-5)*np.ones(len(radii))) * wavelength_U**2 / c / distance**2
    magblackbody = -2.5*np.log10(luminosityblackbody)-48.6
    magblackbody_B = -2.5*np.log10(luminosityblackbody_B)-48.6
    magblackbody_r = -2.5*np.log10(luminosityblackbody_r)-48.6
    magblackbody_i = -2.5*np.log10(luminosityblackbody_i)-48.6
    magblackbody_V = -2.5*np.log10(luminosityblackbody_V)-48.6
    magblackbody_U = -2.5*np.log10(luminosityblackbody_U)-48.6
    # Generate the new curve
    L = ((2.*M*f)/(td)) * (np.exp((-t**2)/td**2)) * E * my_int
    magnitude = -2.5*np.log10(L/4e33)+4.3
    source.data = dict(x=t*(1.+redshift) + T_slider.value, y= magblackbody ,yB = magblackbody_B, yr = magblackbody_r,yi = magblackbody_i, yV = magblackbody_V, yU = magblackbody_U)

for w in [M_slider,f_slider,v_slider, k_slider, T_slider]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = widgetbox(text, M_slider, f_slider, v_slider, k_slider, T_slider)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Sliders"
