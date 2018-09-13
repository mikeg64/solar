
# coding: utf-8

# ## Build a Solar Atmosphere Using McWhirter and VALIIc Data
# 
# ### Steps
# 1. Set parameters for model e.g. physical dimensions and grid size
# 2. Load VALIIIc and McWhirter Data
# 3. Interpolate Data
# 4. Regenerat model
# 

# In[1]:


import os
import numpy as np

import astropy.table
from astropy.table import Table
import astropy.units as u
from astropy.constants import k_B, m_p
import scipy.constants as asc
from astropy.table import Table




#identify location of source data files
#__files__=''
#homedir = os.environ['HOME']
#cwd = os.path.dirname(__files__)
#cwd = homedir+'/Dropbox/multi_tube/python/allpapers/'
#VAL_file = os.path.join(cwd, 'VALIIIC.dat')
#MTW_file = os.path.join(cwd, 'mcwhirter.dat')

#filenames = [VAL_file, MTW_file]



"""
    Read in the data from Table 12 in Vernazza (1981) and combine with
    McWhirter (1975).

    Parameters
    ----------
    VAL_file : string
        The data file for the VAL3c atmosphere, defaults to
        `pysac.mhs_atmosphere.hs_model.VALIIIc_data`

    MTW_file : string
        The data file for the McWhirter atmosphere, defaults to
        `pysac.mhs_atmosphere.hs_model.MTWcorona_data`, if ``False`` is specified
        only the VAL atmosphere is returned.

    mu : float
        The mean molecular weight ratio for the corona. defaults to 0.6.

    Returns
    -------
    data : `astropy.table.Table`
        The combined data, sorted by Z.
"""

def val3c():

	mu=0.602
	VAL_file='VALIIIC.dat'
	VAL3c = Table.read(VAL_file, format='ascii', comment='#',names=['Z','rho','p','T','n_i','n_e','mu'])

	#print VAL3c
	VAL3c['Z'].unit = u.km
	VAL3c['rho'].unit = u.Unit('g cm-3')
	VAL3c['p'].unit = u.Unit('dyne/cm^2')
	VAL3c['T'].unit = u.K
	VAL3c['n_i'].unit = u.one/u.cm**3
	VAL3c['n_e'].unit = u.one/u.cm**3
	# Calculate the mean molecular weight ratio
	VAL3c['mu'] = 4.0/(3*0.74+1+VAL3c['n_e']/VAL3c['n_i'])
	#    VAL3c['mu'] = 4.0/(3*0.74+1+VAL3c['n_e'].quantity/VAL3c['n_i'].quantity)
	#print VAL3c['mu']
	MTW_file='mcwhirter.dat'  
	MTW = Table.read(MTW_file, format='ascii', comment='#', names=['Z', 'T','p','rho'])
	#MTW = Table.read(MTW_file, format='ascii', comment='#')
	#print MTW
	MTW['Z'].unit = u.km
	MTW['T'].unit = u.K
	MTW['p'].unit = u.Unit('dyne cm-2')
	MTW['rho'] = (MTW['p'] / k_B / MTW['T'] * m_p * mu).to('g cm-3')

	#MTW['mu'] = mu

	data = astropy.table.vstack([VAL3c, MTW], join_type='inner')
	#data.sort('Z')
	#    data = astropy.table.vstack([VAL3c, MTW], join_type='inner')
	#print data
	return data   
    

