# import json
# SN2006ca = json.load(open("SN2006ca.json"))
# photometry = SN2006ca["SN2006ca"]["photometry"]

# array = []
# for x in photometry:
# 	for val in x:
# 		if val
# 		print()
		

import requests
import numpy as np
import matplotlib.pyplot as plt

def querry_single_osc(object_name):
    '''
    Download a light curve from the Open Supernova Catalog.
    The output will be in MJD and Magnitude
    '''
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
count = 0
for x in photometry_time:
	#plt.plot(x, photometry_mag[count], "o")
	if photometry_band[count] == "r\'":
		plt.plot(x, photometry_mag[count], "o", color='green')
	elif photometry_band[count] == "i\'":
		plt.plot(x, photometry_mag[count], "o", color='blue')
	elif photometry_band[count] == "V":
		plt.plot(x, photometry_mag[count], "o", color='purple')
	elif photometry_band[count] == "U":
		plt.plot(x, photometry_mag[count], "o", color='orange')
	elif photometry_band[count] == "B":
		plt.plot(x, photometry_mag[count], "o", color='red')
	count += 1

ax = plt.gca()
ax.invert_yaxis()
plt.show()

