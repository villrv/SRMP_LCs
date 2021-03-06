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
from flask import Flask, request, render_template
from bokeh.embed import components

app = Flask(__name__)
@app.route('/')
def index():
    photometry_time, photometry_mag, photometry_sigma, photometry_band, detection = querry_single_osc("DES17C1ffz")

    lumdist, redshift = querry_distance("DES17C1ffz")

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

    callback2=CustomJS(args=dict(source=source), code="""
    	var request = new XMLHttpRequest();
    	request.open('GET', 'https://astrocats.space/api/DES17C1ffz/lumdist+redshift');
    	console.log(request);

    	""")
    callback=CustomJS(args=dict(source=source), code="""

    function numerically_integrate(a, b, dx, f,td) {
        
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
            Ai = dx * (f(x0,td) + f(x1,td))/ 2.;
            
            // cumulatively sum the areas
            Area = Area + Ai    
            
        } 
        return Area;
    }


    //define function to be integrated
    function f1(x,td){
        return x / td * Math.exp(Math.pow(x/td,2)) * Math.exp(-x/8.77);
    }

    function f2(x,td){
        return x / td * Math.exp(Math.pow(x/td,2)) * Math.exp(-x/111);
    }

    dx = 0.5;       // width of the trapezoids

    var tni = 8.8;
    var tco = 111.3;
    var epni = 3.9e10;
    var epco = 6.8e9;
    var beta = 13.8;
    var c = 1.0e10;
    var msol = 2e33;
    var m = mej.value;
    m = m * msol;
    var wav = 6e-5;
    var mni = fni.value;
    mni = mni * m;
    var T = T.value;
    var v = vej.value;
    v = v * Math.pow (10,5);
    var k = k.value;
    var td = Math.sqrt(2. * k * m / (beta * c * v)) / 86400;
    var data = source.data;
    var xstop = 0.0;
    x = data['x']
    y = data['y']
    console.log(x);
    console.log(y);
    var distance = distance.value * 3e24;
    var redshift = redshift.value;
    redshift = parseFloat(redshift);
    for (j = 0; j < x.length; j++) {
        xstop = j * dx;
        int1 = numerically_integrate(0,xstop,dx,f1,td);
        int2 = numerically_integrate(0,xstop,dx,f2,td);
        factor = 2 * mni / td * Math.exp(-Math.pow(xstop/td,2));
        L = factor * ((epni-epco) * int1 + epco * int2) / (4.*3.14*Math.pow(distance,2));
        y[j] = -2.5 * Math.log10(L*wav/c)-48.3;
        x[j] = (dx*(1+redshift)) * j + T;
    }
    source.change.emit();
    """)


    plot = figure(plot_height=400, plot_width=400, title="Super cool blackbody curve thing",
                  tools="crosshair,pan,reset,save,wheel_zoom",
                  x_range=[np.min(photometry_time) - 20, np.max(photometry_time) + 100], y_range=[np.max(photometry_mag), np.min(photometry_mag)])


    plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
    #plot.line('x', 'yB', source=source, line_width=3, line_alpha=0.6, color="pink")
    #plot.line('x', 'yr', source=source, line_width=3, line_alpha=0.6, color="orange")
    #plot.line('x', 'yi', source=source, line_width=3, line_alpha=0.6, color="blue")
    #plot.line('x', 'yV', source=source, line_width=3, line_alpha=0.6, color="turquoise")
    #plot.line('x', 'yU', source=source, line_width=3, line_alpha=0.6, color="purple")

    arrayoftimes = np.array(photometry_time)

    text = TextInput(title="title", value='my parabola', callback = callback2)
    lumdist_input = TextInput(title="title", value=str(lumdist))
    redshift_input = TextInput(title="title", value=str(redshift))


    M_slider = Slider(start=0.1, end=10, value=1, step=.1,
                         title="Ejecta Mass", callback=callback)
    f_slider = Slider(start=0.01, end=1.0, value=0.1, step=.01,
                        title="Nickel Fraction", callback=callback)
    v_slider = Slider(start=5000, end=20000, value=10000, step=1000,
                          title="Ejecta Velocity", callback=callback)
    k_slider = Slider(start=0.1, end=0.4, value=0.2, step=.01,
                           title="Opacity", callback=callback)
    T_slider = Slider(title="Time", value= arrayoftimes.min() - 10,
                start= arrayoftimes.min() - 10, end=arrayoftimes.max() + 10,
                      step= 10,callback=callback)

    callback.args["mej"] = M_slider
    callback.args["fni"] = f_slider
    callback.args["vej"] = v_slider
    callback.args["k"] = k_slider
    callback.args["T"] = T_slider

    callback.args["distance"] = lumdist_input
    callback.args["redshift"] = redshift_input

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
    inputs = widgetbox(text, M_slider, f_slider, v_slider, k_slider, T_slider)
    layout = row(
        plot,
        inputs,
    )

    output_file("NewBokeh.html")

    script,div=components(plot)
    return render_template("blah.html", script=script, div=div)
    show(layout)


if __name__ == '__main__':
    app.run(port=5000, debug = True)

