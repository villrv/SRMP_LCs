import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox, column
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

photometry_time, photometry_mag, photometry_sigma, photometry_band, detection = [0,100],[19,10],[1,1],['r','r'],[]
lumdist, redshift = '',''

N = 500
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

x = t
y = L

filt_data = np.loadtxt('./filters.dat',delimiter=',',dtype='str')
x = t
y = np.zeros((len(L),np.size(filt_data,axis=0)))
column_source = dict(zip(filt_data[:,0], y.T),x=x,y=y)
color_dict = dict(zip(filt_data[:,0], filt_data[:,2]))
for thing in color_dict:
    color_dict[thing] = [color_dict[thing]]

wv_dict = dict(zip(filt_data[:,0], filt_data[:,1]))
for thing in wv_dict:
    wv_dict[thing] = [wv_dict[thing]]

source = ColumnDataSource(data=column_source)
source2 = ColumnDataSource(data=column_source)
color_source = ColumnDataSource(data=color_dict)
wave_source = ColumnDataSource(data=wv_dict)


plot = figure(plot_height=400, plot_width=400, title="Circumstellar Material Model",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[np.min(photometry_time) - 20, np.max(photometry_time) + 100], y_range=[np.max(photometry_mag), np.min(photometry_mag)])


#Hacky solution to get input in
lumdist_input = TextInput(title="", value=str(lumdist))
redshift_input = TextInput(title="", value=str(redshift))


from bokeh.io import output_file, show
from bokeh.models.widgets import Button

javaScript="""
name_vej = vej.value;
name_v = name_v.value;
name_mej = mej.value;
name_k = myk.value;
name_t0 = myt0.value;
name_n = myn.value;
name_rho = myrho.value;
name_r0 = myr0.value;


var today = new Date();
var dd = today.getDate();
var mm = today.getMonth()+1;
var yyyy = today.getFullYear();

if(dd<10) {
    dd = '0'+dd
} 

if(mm<10) {
    mm = '0'+mm
} 

today = mm + '/' + dd + '/' + yyyy;

var name_val = {
        "models": {
            "code":"SNIF",
            "date":today,
            "name":name_v,
            "model":"csm",
            "parameters":{
                        "vejecta":{
                            "latex":"Ejecta Velocity (km/s)",
                            "value":name_vej
                        },
                        "mejecta":{
                            "latex":"Ejecta Mass (Msol)",
                            "value":name_mej
                        },
                        "kappa":{
                            "latex":"Opacity (cm^2/g)",
                            "value":name_k
                        },
                        "texplosion":{
                            "latex":"Explosion Day (MJD)",
                            "value":name_t0
                        },
                        "n":{
                            "latex":"SN Ejecta Power-law Index",
                            "value":name_n
                        },
                         "rho":{
                            "latex":"CSM Inner Density (g / cm^2)",
                            "value":name_rho
                        },
                          "r0":{
                            "latex":"CSM Inner Radius (AU)",
                            "value":name_r0
                        }
                    }
                }
    }

var sampleObject = {};
sampleObject[name_v] = name_val;

const filename = name_v+'.json';
filetext = JSON.stringify(sampleObject, null, "\t");
const blob = new Blob([filetext], { type: 'text/json;charset=utf-8;' });

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename)
} else {
    const link = document.createElement('a')
    link.href = URL.createObjectURL(blob)
    link.download = filename
    link.target = '_blank'
    link.style.visibility = 'hidden'
    link.dispatchEvent(new MouseEvent('click'))
}
"""

button = Button(label="Save Model", button_type="success")



callback=CustomJS(args=dict(source=source,lumdist = lumdist_input,redshift=redshift_input,plotrange = plot.x_range,cd = color_source, wd = wave_source), code="""
t00.start = plotrange.start;
t00.end = plotrange.end;
plotrange.change.emit();

function numerically_integrate(a, b, dx, f,t0,n,s,q,A,betar,betaf,gn) {
 
  // calculate the number of trapezoids
  N = (b - a) / dx;
 
  // define the variable for area
  Area = 0;
 
  //loop to calculate the area of each trapezoid and sum.
  for (i = 1; i <= N; i++) {
     //the x locations of the left and right side of each trapezpoid
     x0 = a + (i-1)*dx;
     x1 = a + i*dx;
    
     // the area of each trapezoid
     Ai = dx * (f(x0,t0,n,s,q,A,betar,betaf,gn) + f(x1,t0,n,s,q,A,betar,betaf,gn))/ 2.;
    
     // cumulatively sum the areas
     Area = Area + Ai  
    
  }
  return Area;
}

function f2(x,t0,n,s,q,A,betar,betaf,gn){
  return 2. * Math.PI * Math.pow(A*gn/q,(5.-n)/(n-s)) * Math.pow(betar,5.-n)*gn * Math.pow((3.-s)/(n-s),3) * Math.pow(x*86400,(2.*n+6.*s-n*s-15.)/(n-s));
}



dx = 0.5;     // width of the trapezoids


var beta = 13.8;
var c = 3.0e10;
var msol = 2e33;
var m = mej.value;
m = m * msol;
var v = vej.value;
v = v * Math.pow (10,5);
var k = k.value;
var data = source.data;
var t00 = t00.value;
var s = 2;
var wav = 6e-5;
var n = nindex.value;
var rho = Math.pow(10,rhoo.value);
var r0 = r0o.value * 1.496e13;
var mcsm = 10.0 * msol;
var esn = 3. * Math.pow(v,2) * m / 10.;
var q = rho * Math.pow(r0,s);
var betaf = 1.0;
var delta = 0;
var betar = 1.0;
var A = 1.0;
var gn = 1.0 / (4.0 * Math.PI * (n - delta)) * Math.pow(2.0 * (5.0 - delta) * (n - 5.0) * esn,(n - 3.) / 2.0) / Math.pow((3.0 - delta) * (n - 3.0) * m,(n - 5.0) / 2.0);  
t0 = k * m / (beta * c * r0);

rcsm = (Math.pow((3.0 - s) /
        (4.0 * Math.PI * q) * mcsm + Math.pow(r0,
            3.0 - s),1.0 / (3.0 - s)))

rph = Math.abs(
       Math.pow(-2.0 * (1.0 - s) /
        (3.0 * k * q) + Math.pow(rcsm,1.0 - s),1.0 /
        (1.0 - s)))

mcsmth = Math.abs(4.0 * Math.PI * q / (3.0 - s) * (
       Math.pow(rph,3.0 - s) - Math.pow(r0, 3.0 - s)))

tfs = (Math.abs((3.0 - s) * Math.pow(q, (3.0 - n) / (n - s)) * Math.pow(A * gn,(s - 3.0) / (n - s)) / (4.0 * Math.PI * Math.pow(betaf,3.0 - s))) ** ( (n - s) / ( (n - 3.0) * (3.0 - s))) * Math.pow(mcsmth,(n - s) / ((n - 3.0) * (3.0 - s))));
trs = Math.pow(v / (betar * Math.pow(A * gn / q,1.0 / (n - s))) * Math.pow(1.0 - (3.0 - n) * m /(4.0 * Math.PI * Math.pow(v,3.0 - n) * gn),1.0 / (3.0 - n)),(n - s) / (s - 3.0));

var xstop = 0.0;
x = data['x'];
y = data['y'];

var distance = parseFloat(lumdist.value) * 3e24;
var redshift = parseFloat(redshift.value);
for (j = 0; j < x.length; j++) {
  xstop = j * dx;
   thing2 = 2. * Math.PI * Math.pow(A*gn/q,(5.-n)/(n-s)) * Math.pow(betar,5.-n)*gn * Math.pow((3.-s)/(n-s),3) * Math.pow(xstop*86400,(2.*n+6.*s-n*s-15.)/(n-s)+1.0)/((2.*n+6.*s-n*s-15.)/(n-s)+1.0)
   factor = 0.5/t0 * Math.exp(-xstop * 86400.0/t0);
   x[j] = (dx*(1+redshift)) * j + t00;
   L = factor * thing2;
   L = L/(4.*Math.PI*Math.pow(distance,2));
   y[j] = -2.5 * Math.log10(L*wav/c)-48.3;
}
for(key in cd.data) {
 y2 = data[key];
   for (j = 0; j < x.length; j++) {
    xstop = j * dx;
  thing2 = 2. * Math.PI * Math.pow(A*gn/q,(5.-n)/(n-s)) * Math.pow(betar,5.-n)*gn * Math.pow((3.-s)/(n-s),3) * Math.pow(xstop*86400,(2.*n+6.*s-n*s-15.)/(n-s)+1.0)/((2.*n+6.*s-n*s-15.)/(n-s)+1.0)
   factor = 0.5/t0 * Math.exp(-xstop * 86400.0/t0);
   L = factor * thing2;
   L = L/(4.*Math.PI*Math.pow(distance,2));
   y2[j] = -2.5 * Math.log10(L*wd.data[key]/c)-48.3;
}
source.change.emit();
}
""")




callback2=CustomJS(args=dict(source=source2,lumdist = lumdist_input,redshift=redshift_input,plotrange = plot.x_range, yplotrange = plot.y_range,cd = color_source, wd = wave_source), code="""
    const url = 'https://astrocats.space/api/' + SNName.value + '/photometry+redshift+lumdist';
    var sourcedata = source.data;
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
    try {
          var data = (myJson[SNName.value]["photometry"]);
    }
    catch(err) {
        window.alert(SNName.value+" doesn't seem to exist in sne.space!");
    }
    var ldist = (myJson[SNName.value]["lumdist"]);
    var z = (myJson[SNName.value]["redshift"]);

    for (key in cd.data){
        var yspec = []
        for (i = 0; i < data.length; i++){
            if ((data[i]['band'] == key) & !('realization' in data[i])){
                yspec.push(parseFloat(data[i]['magnitude']));
                }else{
                yspec.push(NaN);
                }
            }
        sourcedata[key] = yspec;
        source.change.emit();
        }
    
    
    var x = [];
    var y = [];
    for (i = 0; i < data.length; i++){
        if ((!data[i]['upperlimit']) &  !('realization' in data[i])){
        x.push(parseFloat(data[i]['time']));
        y.push(parseFloat(data[i]['magnitude']));
        }
    }

    function bouncer(arr){
    return arr.filter(Boolean);
    }

    x = bouncer(x);
    y = bouncer(y);

    sourcedata["x"] = x;
    sourcedata["y"] = y;
    console.log(lumdist.value)
    lumdist.value = ldist[0]['value'];
    redshift.value = z[0]['value'];

    plotrange.start = Math.min(x);
    plotrange.start = Math.min.apply(Math, x);
    yplotrange.start = Math.max.apply(Math, y);
    plotrange.end  = Math.max.apply(Math, x);
    yplotrange.end = Math.min.apply(Math, y);
    plotrange.change.emit();
    yplotrange.change.emit();
    source.change.emit();
    lumdist.change.emit();
    redshift.change.emit();

  })
    """)


plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
for i,filt in enumerate(filt_data[:,0]):
    plot.line(x='x', y=filt, source=source, line_width=3, line_alpha=0.6, 
        color=filt_data[i,2])
plot.circle('x', 'y', source=source2)
for i,filt in enumerate(filt_data[:,0]):
    plot.circle(x='x', y=filt, source=source2, line_width=3, line_alpha=0.6, 
        color=filt_data[i,2])


arrayoftimes = np.array(photometry_time)

text = TextInput(title="Insert the name of the supernova here:", value='', callback = callback2)

M_slider = Slider(start=0.1, end=10, value=1, step=.1,
                     title="Ejecta Mass (solar mass)", callback=callback)
n_slider = Slider(start=6, end=14, value=7, step=0.1,
                    title="Ejecta Shape n", callback=callback)
rho_slider = Slider(start=-15, end=-1, value=-10, step=.1,
                    title="log Density (g/cm^3)", callback=callback)
r0_slider = Slider(start=1, end=100, value=5, step=1,
                    title="CSM Radius (AU)", callback=callback)
v_slider = Slider(start=5000, end=20000, value=10000, step=1000,
                      title="Ejecta Velocity (km / s)", callback=callback)
k_slider = Slider(start=0.1, end=0.4, value=0.2, step=.01,
                       title="Opacity (cm^2 / g)", callback=callback)
T_slider = Slider(title="Time (MJD)", value= arrayoftimes.min() - 10,
            start= arrayoftimes.min() - 10, end=arrayoftimes.max() + 10,
                  step= 10,callback=callback)

callback.args["mej"] = M_slider
callback.args["vej"] = v_slider
callback.args["k"] = k_slider
callback.args["nindex"] = n_slider
callback.args["rhoo"] = rho_slider
callback.args["r0o"] = r0_slider
callback.args["t00"] = T_slider

button.callback = CustomJS(args=dict(source=source,vej=v_slider,
                name_v=text,mej = M_slider, myk = k_slider, 
                myt0 = T_slider, myn = n_slider,
                myrho = rho_slider, myr0 = r0_slider),
                code=javaScript)


callback2.args["SNName"] = text


text.on_change('value', update_title)


# Set up layouts and add to document
inputs = widgetbox(text, M_slider, v_slider, k_slider, n_slider, r0_slider, rho_slider, T_slider)
layout = row(
    plot,
    column(
    inputs,
    button)
)

output_file("CSM.html")

save(plot)
show(layout)

