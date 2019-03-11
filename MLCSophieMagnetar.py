import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox
from bokeh.models import ColumnDataSource, CustomJS
from bokeh.models.widgets import Slider, TextInput
from bokeh.plotting import figure, output_file, save
from scipy import integrate
import requests
from bokeh.plotting import figure, output_file, show
from QuerySingleOSC import querry_single_osc
from UpdateTitle import update_title
from BlackbodyFunction import blackbody
from QueryDistance import querry_distance
from bokeh.models.sources import AjaxDataSource

photometry_time, photometry_mag, photometry_sigma, photometry_band, detection = querry_single_osc("DES17C1ffz")

lumdist, redshift = querry_distance("DES17C1ffz")

N = 200
t_original = np.linspace(0, 500, N)
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


source = ColumnDataSource(data=dict(x=t, y=L, yB=L, yR=L, yi=L, yG=L))


x=t
y=L

source2 = ColumnDataSource(data=dict(x=x, y=y, yG=y, yR=y, yi=y, yB=y))


plot = figure(plot_height=400, plot_width=400, title="Super cool blackbody curve thing",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[np.min(photometry_time) - 20, np.max(photometry_time) + 100], y_range=[np.max(photometry_mag), np.min(photometry_mag)])

callback=CustomJS(args=dict(source=source, plotrange=plot.x_range), code="""
T.start = plotrange.start;
T.end = plotrange.end;
plotrange.change.emit();

function numerically_integrate(a, b, dx, f,td,yy) {
   
   // calculate the number of trapezoids
   n = (b - a) / dx;
   
   // define the variable for area
   Area = 0;
   
   //loop to calculate the area of each trapezoid and sum.
   for (i = 1; i <= n; i++) {
      //the x locations of the left and right side of each trapezpoid
      x0 = a + (i-1)*dx;
      x1 = a + i*dx;
      
      // the area of each trapezoid
      Ai = dx * (f(x0,td,yy) + f(x1,td,yy))/ 2.;
      
      // cumulatively sum the areas
      Area = Area + Ai   
      
   } 
   return Area;
}


//define function to be integrated
function f1(x,td,yy){
   return Math.exp(Math.pow(x/td,2)) * x/td * 1./Math.pow(1. + yy * x/td,2);
}

dx = 0.5;     // width of the trapezoids

var T = T.value;
var beta = 13.8;
var c = 1.0e10;
var msol = 2e33;
var m = mej.value;
m = m * msol;
var v = vej.value;
v = v * Math.pow (10,5);
var wav = 6e-5;
var k = k.value;
var b = bfield.value;
var p = pspin.value;
var td = Math.sqrt(2. * k * m / (beta * c * v)) / 86400;
var data = source.data;
x = data['x']
y = data['y']
var ep = Math.pow(0.1 * p,-2) * (2e50);
var tp = Math.pow(b,-2) * 1.3 * Math.pow(p*0.1,2) * 365.0;
var yy = td/tp;
var distance = distance.value * 3e24;
var redshift = redshift.value;


for (j = 0; j < x.length; j++) {
    xstop = j * dx;
    int1 = numerically_integrate(0,xstop,dx,f1,td,yy);
    factor = 2. * ep / (tp * 86400.) * Math.exp(-Math.pow(xstop/td,2));
    L = factor * int1 / (4.*3.14*Math.pow(distance,2));
    y[j] = -2.5 * Math.log10(L*wav/c)-48.3;
    x[j] = (dx*(1+redshift)) * j + T;
    console.log(int1);
}
source.change.emit();
""")


#l input L_input = (Ep/Tp) * (1/(1+t/Tp)**2)

#e knot m*(v**2)/4 td = (2*k*m/beta*c*v)**(1/2)

callback2=CustomJS(args=dict(source=source2, plotrange=plot.x_range, plotrange_y=plot.y_range), code="""
	var sourcedata = source.data;
	const url = 'https://astrocats.space/api/' + supernovae.value + '/photometry/time+magnitude+e_magnitude+band+upperlimit';
fetch(url, {
  method: "GET",
  mode: 'cors',
  headers: {
    "Content-Type": "text/html"
  }
}).then(function(response) {
    return response.json();
  })
  .then(function(myJson) {
    var data = myJson[supernovae.value]["photometry"];
    var x = [];
    var y = [];
    var yG = [];
    var yR = [];
    var yi = [];
    var yB =[];
    for (i = 0; i < data.length; i++) {
    	if (!data[i][4]) {
    	x.push(parseFloat(data[i][0]));
    	y.push(parseFloat(data[i][1]));
    	if (data[i][3] == "g") {
    	yG.push(parseFloat(data[i][1]));
    	}
    	else {
    	yG.push(NaN);
    	}
    	if (data[i][3] == "r") {
    	yR.push(parseFloat(data[i][1]));
    	}
    	else {
    	yR.push(NaN);
    	}
    	if (data[i][3] == "i") {
    	yi.push(parseFloat(data[i][1]));
    	}
    	else {
    	yi.push(NaN);
    	}
    	if (data[i][3] == "B") {
    	yB.push(parseFloat(data[i][1]));
    	}
    	else {
    	yB.push(NaN);
    	}
    	}

    } 
    function bouncer(arr) {
    return arr.filter(Boolean);
    }
    x = bouncer(x);
    y = bouncer(y);
    sourcedata["x"] = x;
    sourcedata["y"] = y;
    sourcedata["yG"] = yG;
    sourcedata["yR"] = yR;
    sourcedata["yi"] = yi;
    sourcedata["yB"] = yB;
    plotrange.start = Math.min.apply(Math, x);
    plotrange_y.start = Math.max.apply(Math, y);
    plotrange.end  = Math.max.apply(Math, x);
    plotrange_y.end = Math.min.apply(Math, y);
    plotrange.change.emit();
    plotrange_y.change.emit();
    source.change.emit();
  })
	""")


plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6, color="black")
plot.line('x', 'yG', source=source, line_width=3, line_alpha=0.6, color="green")
plot.line('x', 'yR', source=source, line_width=3, line_alpha=0.6, color="red")
plot.line('x', 'yi', source=source, line_width=3, line_alpha=0.6, color="purple")
plot.line('x', 'yB', source=source, line_width=3, line_alpha=0.6, color="blue")
plot.circle('x', 'y', source=source2, size = 5, color="black")
plot.circle('x', 'yG', source=source2, size = 5, color="green")
plot.circle('x', 'yB', source=source2, size = 5, color="blue")
plot.circle('x', 'yR', source=source2, size = 5, color="red")
plot.circle('x', 'yi', source=source2, size = 5, color="purple")
#plot.line('x', 'yV', source=source, line_width=3, line_alpha=0.6, color="turquoise")
#plot.line('x', 'yU', source=source, line_width=3, line_alpha=0.6, color="purple")

arrayoftimes = np.array(photometry_time)

text = TextInput(title="title", value='my parabola', callback = callback2)
lumdist_input = TextInput(title="title", value=str(lumdist))
redshift_input = TextInput(title="title", value=str(redshift))


M_slider = Slider(start=0.1, end=10, value=1, step=.1,
                     title="Ejecta Mass", callback=callback)
v_slider = Slider(start=5000, end=20000, value=10000, step=1000,
                      title="Ejecta Velocity", callback=callback)
k_slider = Slider(start=0.1, end=0.4, value=0.2, step=.01,
                       title="Opacity", callback=callback)
T_slider = Slider(title="Time", value= arrayoftimes.min() - 10,
            start= arrayoftimes.min() - 10, end=arrayoftimes.max() + 10,
                  step= 10,callback=callback)
bfield_slider = Slider(start=1.0, end=10.0, value=5, step=1.0,
                    title="Ep", callback=callback)
pspin_slider = Slider(start=1.0, end=10.0, value=5, step=1.0,
                    title="Tp", callback=callback)

callback.args["mej"] = M_slider
callback.args["vej"] = v_slider
callback.args["k"] = k_slider
callback.args["T"] = T_slider
callback.args["bfield"] = bfield_slider
callback.args["pspin"] = pspin_slider
callback.args["distance"] = lumdist_input
callback.args["redshift"] = redshift_input



callback2.args["supernovae"] = text


plot.x_range.start = 10

# count = 0
# for x in photometry_time:
# 	if photometry_band[count] == "r":
# 		plot.circle(x, photometry_mag[count], size=5, color="orange", alpha=0.5)
# 	elif photometry_band[count] == "i":
# 		plot.circle(x, photometry_mag[count], size=5, color="blue", alpha=0.5)
# 	elif photometry_band[count] == "g":
# 		plot.circle(x, photometry_mag[count], size=5, color="turquoise", alpha=0.5)
# 	elif photometry_band[count] == "z":
# 		plot.circle(x, photometry_mag[count], size=5, color="purple", alpha=0.5)
# 	elif photometry_band[count] == "B":
# 		plot.circle(x, photometry_mag[count], size=5, color="pink", alpha=0.5)
# 	count += 1

text.on_change('value', update_title)

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


#for w in [M_slider,f_slider,v_slider, k_slider, T_slider]:
#    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = widgetbox(text, M_slider, v_slider, k_slider, T_slider, bfield_slider, pspin_slider)
layout = row(
    plot,
    inputs,
)

output_file("NewBokehMagnetar.html")

save(plot)
show(layout)

