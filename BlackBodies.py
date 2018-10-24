import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput
from bokeh.plotting import figure
from scipy import integrate
import math

# Set up data
sigb = 5.67*10**-8
c = 299792458
h = 6.62*10**-34
k = 1.3*10**-23
wav = 5e-5
T = 5000
firstpart = (2*math.pi*c**2*h) / wav
complexthing = math.exp**(((h*c)/(wav*k*T)) - 1)
secondpart = 1 / complexthing
F = firstpart*secondpart

source = ColumnDataSource(data=dict(x=T, y=F))

# Set up plot
plot = figure(plot_height=400, plot_width=400, title="Black Body",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0, 100], y_range=[-15, -19])

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)

def update_title(attrname, old, new):
    plot.title.text = text.value

text.on_change('value', update_title)

source.data = dict(x=T, y=F)

curdoc().add_root(row(inputs, plot, width=800))