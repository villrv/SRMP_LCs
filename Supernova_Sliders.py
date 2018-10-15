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
# x = np.array([0, 1, 2, 3, 4, 5])
td = 10
m = 5
f = 0.5
integrand = x/td*(np.exp(x**2/td**2))*(np.exp(-x/td))
y_int = integrate.cumtrapz(integrand, x, initial = 0)
y = y_int*2*m*f/td
source = ColumnDataSource(data=dict(x=x, y=y))


# Set up plot
plot = figure(plot_height=400, plot_width=400, title="Cubic",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0, 100], y_range=[0, 1000])

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)


# Set up widgets
text = TextInput(title="Title", value='My Cubic Function')
massejecta = Slider(title="Mass of ejecta", value=1.0, start=0.0, end=10.0, step=1)
fracradioactive = Slider(title="Radioactive stuff", value=0.5, start=0.0, end=1.0, step=0.1)
diffusiontime = Slider(title="Diffusion time", value=10, start=1.0, end=100.0, step=1)

WeirdThingyThatFlattensTheCurve = Slider(title="Weird thingy that flattens the curve", value=0.0, start=0.0, end=2*np.pi)
height = Slider(title="Height", value=1.0, start=0.1, end=5.1, step=0.1)
# integrand = x**2


# Set up callbacks
def update_title(attrname, old, new):
    plot.title.text = text.value

text.on_change('value', update_title)

def update_data(attrname, old, new):

    # Get the current slider values
    f = fracradioactive.value
    m = massejecta.value
    td = diffusiontime.value

    w = WeirdThingyThatFlattensTheCurve.value
    k = height.value

    # myint = cumtrapz(x, integrand)

    # Generate the new curve
    x = np.linspace(0, 100, N)
    integrand = x/td*(np.exp(x**2/td**2))*(np.exp(-x/td))
    y_int = integrate.cumtrapz(integrand, x, initial = 0)
    y = y_int*2*m*f/td
    print(y)

    source.data = dict(x=x, y=y)

for w in [fracradioactive, massejecta, diffusiontime, WeirdThingyThatFlattensTheCurve, height]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = widgetbox(text, fracradioactive, massejecta, diffusiontime, WeirdThingyThatFlattensTheCurve, height)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Sliders"
