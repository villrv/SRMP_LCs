''' Present an interactive function explorer with slider widgets.

Scrub the sliders to change the properties of the ``sin`` curve, or
type into the title text box to update the title of the plot.

Use the ``bokeh serve`` command to run the example by executing:

    bokeh serve sliders.py

at your command prompt. Then navigate to the URL

    http://localhost:5006/sliders

in your browser.

'''
import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput
from bokeh.plotting import figure
from scipy import integrate
# From SN2006 v
import json
sn = json.load(open("SN2006ca.json"))
photom = print(sn["SN2006ca"]["photometry"])
import requests
import math

# Set up data
N = 200
x = np.linspace(0, 100, N)
td = 10
m = 5
f = 0.5
sig = 5.67e-5
eval = (np.exp(x**2/td**2))
integrand = x/td*eval*(np.exp(-x/td))
y_int = integrate.cumtrapz(integrand, x, initial = 0)
y = y_int*2*m*f/td
source = ColumnDataSource(data=dict(x=x, y=y, yU=y, yV=y, yB=y, yR=y, yI=y))



# Set up plot ; also 4,320,000 seconds is 50 days
plot = figure(plot_height=400, plot_width=400, title="Supernova",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0, 100], y_range=[-19, 19])

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)

plot.line('x', 'yU', source=source, line_width=3, line_alpha=0.6)

plot.line('x', 'yV', source=source, line_width=3, line_alpha=0.6)

plot.line('x', 'yB', source=source, line_width=3, line_alpha=0.6)

plot.line('x', 'yR', source=source, line_width=3, line_alpha=0.6)

plot.line('x', 'yI', source=source, line_width=3, line_alpha=0.6)


# Set up widgets ; replace td with velocity and opacity
text = TextInput(title="Title", value='My Supernova')
massejecta = Slider(title="Mass of ejecta (super big)", value=5, start=1, end=20, step=1)
fracradioactive = Slider(title="Radioactive stuff", value=0.5, start=0.0, end=1.0, step=0.1)
# diffusiontime = Slider(title="Diffusion time", value=10, start=1.0, end=200.0, step=1)
velocity = Slider(title="Velocity", value=10000, start=5000, end=20000, step=1000)
opacity = Slider(title="Opacity", value=0.2, start=0.1, end=0.5, step=0.1)

# WeirdThingyThatFlattensTheCurve = Slider(title="Weird thingy that flattens the curve", value=0.0, start=0.0, end=2*np.pi)
# height = Slider(title="Height", value=1.0, start=0.1, end=5.1, step=0.1)
# integrand = x**2


# Set up callbacks
def update_title(attrname, old, new):
    plot.title.text = text.value

text.on_change('value', update_title)

def blackbody(r, T, wav): 
    sigb = 5.67*10**-8
    c = 3e10
    h = 6.26e-27
    k = 1.38e-16
    firstpart = (2*math.pi*c**2*h) / wav**5
    complexthing = np.exp(((h*c)/(wav*k*T)) - 1.)
    secondpart = 1. / complexthing
    F = firstpart*secondpart
    flux = F*r**2
    return flux

def update_data(attrname, old, new):

    # Get the current slider values
    f = fracradioactive.value
    m = massejecta.value * 2.e33
    v = velocity.value * 1.e5
    k = opacity.value
    e = 3.9*10**10 # epsilon ... some constant
    c = 3.*10**10 # c is the speed of light
    b = 13.8
    tn = 8.8
    td = (2.*k*m/(b*c*v))**0.5/86400.


    # Generate the new curve
    x = np.linspace(0, 100, N)
    integrand = x/td*(np.exp(x**2/td**2))*(np.exp(-x/tn))
    y_int = integrate.cumtrapz(integrand, x, initial = 0)
    neg = (np.exp(-x**2/td**2))
    y = e*neg*2*m*f/(td)*y_int
    mag = -2.5*np.log10(y/4e33) + 4.3
    rad = v*x*86400
    temp = (y/(4.*math.pi*sig*rad**2))**0.25
    wav = 5e-5
    wavU = 3.65e-5
    wavV = 5.51e-5
    wavB = 4.45e-5
    wavR = 6.58e-5
    wavI = 8.06e-5
    d = 8e26
    lboriginal = blackbody(rad, temp, wav*np.ones(len(rad)))
    lboriginal = lboriginal*wav**2/(c*d**2)
    magbb = -2.5*np.log10(lboriginal)-48.6

    lbU = blackbody(rad, temp, wavU*np.ones(len(rad)))
    lbU = lbU*wavU**2/(c*d**2)
    magU = -2.5*np.log10(lbU)-48.6

    lbV = blackbody(rad, temp, wavV*np.ones(len(rad)))
    lbV = lbV*wavV**2/(c*d**2)
    magV = -2.5*np.log10(lbV)-48.6

    lbB = blackbody(rad, temp, wavB*np.ones(len(rad)))
    lbB = lbB*wavB**2/(c*d**2)
    magB = -2.5*np.log10(lbB)-48.6

    lbR = blackbody(rad, temp, wavR*np.ones(len(rad)))
    lbR = lbR*wavR**2/(c*d**2)
    magR = -2.5*np.log10(lbR)-48.6

    lbI = blackbody(rad, temp, wavI*np.ones(len(rad)))
    lbI = lbI*wavI**2/(c*d**2)
    magI = -2.5*np.log10(lbI)-48.6

    print(magbb)


    # print(mag)

    source.data = dict(x=x, y=magbb, yU=magU, yV=magV, yB=magB, yR=magR, yI=magI)

    # source.data = dict(x=x, yU=magU)

    # source.data = dict(x=x, yV=magV)

    # source.data = dict(x=x, yB=magB)

    # source.data = dict(x=x, yR=magR)

    # source.data = dict(x=x, yI=magI)

for w in [fracradioactive, massejecta, velocity, opacity]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = widgetbox(text, fracradioactive, massejecta, velocity, opacity)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Sliders"

# !! From SN2006 !!

# plt.style.use('seaborn-whitegrid')

def querry_single_osc(object_name):
    
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

x = np.zeros(len(photometry_time))
y = np.zeros(len(photometry_mag))

x = photometry_time
x = x - 53800
y = photometry_mag

s = photometry_sigma
b = photometry_band
# l = photometry_limit

for i in np.arange(len(x)):
    if photometry_band[i] == "U":
        plot.circle([x[i]], [y[i]], size=20, color="purple", alpha=0.5)
    elif photometry_band[i] == "B":
        plot.circle([x[i]], [y[i]], size=20, color="pink", alpha=0.5)
    elif photometry_band[i] == "V":
        plot.circle([x[i]], [y[i]], size=20, color="blue", alpha=0.5)
    elif photometry_band[i] == "r":
        plot.circle([x[i]], [y[i]], size=20, color="yellow", alpha=0.5)
    elif photometry_band[i] == "i":
        plot.circle([x[i]], [y[i]], size=20, color="black", alpha=0.5)

# print(x)

# plt.gca().invert_yaxis()
# plt.show()
# show(plot)
