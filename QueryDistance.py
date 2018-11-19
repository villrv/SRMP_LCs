import numpy as np
import requests

def querry_distance(object_name):
    '''
    Download a light curve from the Open Supernova Catalog.
    The output will be in MJD and Magnitude
    '''
    osc_link    = 'https://astrocats.space/api/' + object_name + '/lumdist+redshift'
    osc_request = requests.get(osc_link).json()
    osc_data    = osc_request[object_name]
    lumdist = np.array(osc_data['lumdist'][0]['value']).astype(float)
    redshift = np.array(osc_data['redshift'][0]['value']).astype(float)

    return lumdist, redshift