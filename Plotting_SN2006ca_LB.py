import json
sn = json.load(open("SN2006ca.json"))
photom = print(sn["SN2006ca"]["photometry"])
import requests
import numpy as np
import matplotlib.pyplot as plt

# matplotlib inline
plt.style.use('seaborn-whitegrid')

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
y = photometry_mag

s = photometry_sigma
b = photometry_band
# l = photometry_limit

for i in np.arange(len(x)):
	if photometry_band[i] == "U":
		plt.plot(x[i], y[i], 'o', color = 'purple')
	elif photometry_band[i] == "B":
		plt.plot(x[i], y[i], 'o', color = 'pink')
	elif photometry_band[i] == "V":
		plt.plot(x[i], y[i], 'o', color = 'blue')
	elif photometry_band[i] == "r":
		plt.plot(x[i], y[i], 'o', color = 'yellow')
	elif photometry_band[i] == "i":
		plt.plot(x[i], y[i], 'o', color = 'black')

plt.gca().invert_yaxis()
plt.show()
