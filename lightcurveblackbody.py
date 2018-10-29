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

photometry_time, photometry_mag, photometry_sigma, photometry_band, detection = querry_single_osc("SN2006ca")


N = 200
t = np.linspace(0, 100, N)
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


source = ColumnDataSource(data=dict(x=t, y=L))


# Set up plot
plot = figure(plot_height=400, plot_width=400, title="super cool parabola",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0, 100], y_range=[10**44, 2*10**45])

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)


# Set up widgets
text = TextInput(title="title", value='my parabola')
M_slider = Slider(title="Ejecta Mass", value=5, start=1.0, end=10.0, step=0.5)
f_slider = Slider(title="Fraction of Radio Active Stuff", value=0.0, start=0.0, end=1.0, step=0.1)
v_slider = Slider(title="Ejecta Velocity", value= 12000, start=5000, end=20000, step=1000)
k_slider = Slider(title="Kapa", value=0.25, start= 0.1, end=0.5, step=0.05)

count = 0
for x in photometry_time:
	x = x - 53860
	#plt.plot(x, photometry_mag[count], "o")
	if photometry_band[count] == "r\'":
		plot.circle(x, photometry_mag[count], size=5, color="green", alpha=0.5)
	elif photometry_band[count] == "i\'":
		plot.circle(x, photometry_mag[count], size=5, color="blue", alpha=0.5)
	elif photometry_band[count] == "V":
		plot.circle(x, photometry_mag[count], size=5, color="purple", alpha=0.5)
	elif photometry_band[count] == "U":
		plot.circle(x, photometry_mag[count], size=5, color="orange", alpha=0.5)
	elif photometry_band[count] == "B":
		plot.circle(x, photometry_mag[count], size=5, color="red", alpha=0.5)
	count += 1

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
    wavelength = 5*10**-5
    distance = 3.*(10**19)
    luminosityblackbody = blackbody(radii, temperature, 5.*(10**-5)*np.ones(len(radii))) * wavelength**2 / c / distance**2
    magblackbody = -2.5*np.log10(luminosityblackbody)-48.6
    print(magblackbody)
    # Generate the new curve
    L = ((2.*M*f)/(td)) * (np.exp((-t**2)/td**2)) * E * my_int
    magnitude = -2.5*np.log10(L/4e33)+4.3

    source.data = dict(x=t, y=luminosityblackbody)

for w in [M_slider,f_slider,v_slider, k_slider]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = widgetbox(text, M_slider, f_slider, v_slider, k_slider)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Sliders"
