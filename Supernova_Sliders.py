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

# Set up data
N = 200
x = np.linspace(0, 100, N)
x = x
# x = np.array([0, 1, 2, 3, 4, 5])
td = 10
m = 5
f = 0.5
eval = (np.exp(x**2/td**2))
integrand = x/td*eval*(np.exp(-x/td))
y_int = integrate.cumtrapz(integrand, x, initial = 0)
y = y_int*2*m*f/td
source = ColumnDataSource(data=dict(x=x, y=y))


# Set up plot ; also 4,320,000 seconds is 50 days
plot = figure(plot_height=400, plot_width=400, title="Supernova",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0, 100], y_range=[0, 10**43])

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)


# Set up widgets ; replace td with velocity and opacity
text = TextInput(title="Title", value='My Supernova')
massejecta = Slider(title="Mass of ejecta (super big)", value=5*10**33, start=2*10**33, end=20*10**33, step=1*10**33)
fracradioactive = Slider(title="Radioactive stuff", value=0.5, start=0.0, end=1.0, step=0.1)
# diffusiontime = Slider(title="Diffusion time", value=10, start=1.0, end=200.0, step=1)
velocity = Slider(title="Velocity", value=10*10**8, start=5*10**8, end=20*10**8, step=1*10**8)
opacity = Slider(title="Opacity", value=0.2, start=0.1, end=0.5, step=0.1)

# WeirdThingyThatFlattensTheCurve = Slider(title="Weird thingy that flattens the curve", value=0.0, start=0.0, end=2*np.pi)
# height = Slider(title="Height", value=1.0, start=0.1, end=5.1, step=0.1)
# integrand = x**2


# Set up callbacks
def update_title(attrname, old, new):
    plot.title.text = text.value

text.on_change('value', update_title)

def update_data(attrname, old, new):

    # Get the current slider values
    f = fracradioactive.value
    m = massejecta.value
    v = velocity.value
    k = opacity.value
    e = 3.9*10**10 # epsilon ... some constant
    c = 3*10**10 # c is the speed of light
    b = 13.8
    tn = 8.8
    td = (2*k*m/(b*c*v))**0.5/84600.



    # Generate the new curve
    x = np.linspace(0, 100, N)
    integrand = x/td*(np.exp(x**2/td**2))*(np.exp(-x/tn))
    y_int = integrate.cumtrapz(integrand, x, initial = 0)
    neg = (np.exp(-x**2/td**2))
    y = e*neg*2*m*f/(86400.*td)*y_int
    print(y)

    source.data = dict(x=x, y=y)

for w in [fracradioactive, massejecta, velocity, opacity]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = widgetbox(text, fracradioactive, massejecta, velocity, opacity)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Sliders"
