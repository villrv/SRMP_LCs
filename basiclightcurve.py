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

from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput
from bokeh.plotting import figure

# Set up data
N = 200
t = np.linspace(0, 100, N)
M = 1
f = 0.5
td = 10
E = 1.0
integrand = t**2
my_int = integrate.cumtrapz(t,integrand, initial = 0)
L = ((2*M*f)/td) * 2.71828**((-t**2)/(td**2)) * E * my_int
print(L)

source = ColumnDataSource(data=dict(x=t, y=L))


# Set up plot
plot = figure(plot_height=400, plot_width=400, title="super cool parabola",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0, 100], y_range=[0, 100])

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)


# Set up widgets
text = TextInput(title="title", value='my parabola')
M_slider = Slider(title="a value", value=0.0, start=-0.0, end=10.0, step=1.0)
f_slider = Slider(title="b value", value=0.0, start=-0.0, end=1.0, step=0.1)
td_slider = Slider(title="c value", value=1.0, start=-1.0, end=100.0, step=5.0)


# Set up callbacks
def update_title(attrname, old, new):
    plot.title.text = text.value

text.on_change('value', update_title)

def update_data(attrname, old, new):

    # Get the current slider values
    Mvalue = M_slider.value
    fvalue = f_slider.value
    tdvalue = td_slider.value

    # Generate the new curve
    L = ((2*M*f)/td) * 2.71828**((-t**2)/(td**2)) * E * my_int

    source.data = dict(x=t, y=L)

for tdvalue in [M_slider,f_slider,td_slider]:
    tdvalue.on_change('value', update_data)


# Set up layouts and add to document
inputs = widgetbox(text, M_slider, f_slider, td_slider)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Sliders"
