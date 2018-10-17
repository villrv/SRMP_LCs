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
              x_range=[0, 4320000], y_range=[0, 10**44])

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)


# Set up widgets
text = TextInput(title="title", value='my parabola')
M_slider = Slider(title="Ejecta Mass", value=5*(1.989*(10**33)), start=1.0*(1.989*(10**33)), end=10.0*(1.989*(10**33)), step=0.5*(1.989*(10**33)))
f_slider = Slider(title="Fraction of Radio Active Stuff", value=0.0, start=0.0, end=1.0, step=0.1)
v_slider = Slider(title="Ejecta Velocity", value= 12*(10**8), start=5*(10**8), end=20*(10**8), step=1*(10**8))
k_slider = Slider(title="Kapa", value=0.25, start= 0.1, end=0.5, step=0.05)


# Set up callbacks
def update_title(attrname, old, new):
    plot.title.text = text.value

text.on_change('value', update_title)

def update_data(attrname, old, new):

    # Get the current slider values
    M = M_slider.value
    f = f_slider.value
    v = v_slider.value
    k = k_slider.value
    tn = 8.8
    B = 13.8
    c = 3*(10**10)
    E = 3.9*(10**10)
    td = (((2*k*M)/(B*c*v))**(1./2.))/60./60./24.
    integrand = (t/td)*(np.exp(t**2/td**2))*(np.exp(-t/tn))
    my_int = integrate.cumtrapz(t,integrand, initial = 0)

    # Generate the new curve
    L = ((2*M*f)/(td*24.*60.*60.)) * (np.exp((-t**2)/td**2)) * E * my_int
    print(L)

    source.data = dict(x=t, y=L)

for w in [M_slider,f_slider,v_slider, k_slider]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = widgetbox(text, M_slider, f_slider, v_slider, k_slider)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Sliders"
