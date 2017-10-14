# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 11:37:39 2014

@author: sm1fg

    Generate a 1D non-magnetic atmosphere vector based on an empirical model 
    based on observational data, or specify an analytical hydrostatic 
    equilibrium atmosphere.

"""
import os
import numpy as np
import astropy.units as u
from scipy.interpolate import UnivariateSpline


# Fibonacci numbers module

def fib(n):    # write Fibonacci series up to n
    a, b = 0, 1
    while b < n:
        print(b, end=' ')
        a, b = b, a+b
    print()

def fib2(n):   # return Fibonacci series up to n
    result = []
    a, b = 0, 1
    while b < n:
        result.append(b)
        a, b = b, a+b
    return result




#============================================================================
# Read in and interpolate HD atmosphere
#============================================================================

def read_VAL3c_MTW(VAL_file=None, MTW_file=None, mu=0.602):
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
    from . import VALIIIc_data, MTWcorona_data
    if not VAL_file:
        VAL_file = VALIIIc_data
    if MTW_file is None:
        MTW_file = MTWcorona_data

    VAL3c = Table.read(VAL_file, format='ascii', comment='#')
    VAL3c['Z'].unit = u.km
    VAL3c['rho'].unit = u.Unit('g cm-3')
    VAL3c['p'].unit = u.Unit('dyne/cm^2')
    VAL3c['T'].unit = u.K
    VAL3c['n_i'].unit = u.one/u.cm**3
    VAL3c['n_e'].unit = u.one/u.cm**3
    # Calculate the mean molecular weight ratio
    VAL3c['mu'] = 4.0/(3*0.74+1+VAL3c['n_e']/VAL3c['n_i'])
#    VAL3c['mu'] = 4.0/(3*0.74+1+VAL3c['n_e'].quantity/VAL3c['n_i'].quantity)

    if MTW_file:
        MTW = Table.read(MTW_file, format='ascii', comment='#')
        MTW['Z'].unit = u.km
        MTW['T'].unit = u.K
        MTW['p'].unit = u.Unit('dyne cm-2')
        MTW['rho'] = (MTW['p'] / k_B / MTW['T'] * m_p * mu).to('g cm-3')

        MTW['mu'] = mu

        data = astropy.table.vstack([VAL3c, MTW], join_type='inner')
        #    data = astropy.table.vstack([VAL3c, MTW], join_type='inner')

    else:
        data = VAL3c

    data.sort('Z')

    return data

def read_dalsgaard(DAL_file=None, mu=0.602):
    """
    Read in the data from Table in Christensen-Dalsgaard (1996).

    Parameters
    ----------
    DAL_file : string
        The data file for the VAL3c atmosphere, defaults to
        `pysac.mhs_atmosphere.hs_model.VALIIIc_data`

    mu : float
        The mean molecular weight ratio for solar interior. defaults to 0.602
        for fully ionized plasma.

    Returns
    -------
    data : `astropy.table.Table`
        The combined data, sorted by Z.
    """
    from . import dalsgaard_data
    if not DAL_file:
        DAL_file = dalsgaard_data

    DAL = Table.read(DAL_file, format='ascii', comment='#')
    DAL['Z'] *= 6.96342e8 # convert from ratio of solar radius to m
    DAL['Z'].unit = u.m
    DAL['sound_speed'].unit = u.Unit('cm/s')
    DAL['rho'].unit = u.Unit('g cm-3')
    DAL['p'].unit = u.Unit('dyne/cm^2')
    DAL['T'].unit = u.K
    DAL['Gamma_1'].unit = u.one

    # Calculate the mean molecular weight ratio
    #VAL3c['mu'] = 4.0/(3*0.74+1+VAL3c['n_e']/VAL3c['n_i'])

    data = astropy.table.vstack([VAL3c, MTW], join_type='inner')
#    data = astropy.table.vstack([VAL3c, MTW], join_type='inner')
    data.sort('Z')

    return data

#============================================================================
# interpolate the empirical data onto a Z array
#============================================================================
def interpolate_atmosphere(data, Z, s=0.25):
    """ This module generates a 1d array for the model plasma preesure, plasma
    density, temperature and mean molecular weight.
    """

    from scipy.interpolate import UnivariateSpline
    hdata = np.array(u.Quantity(data['Z']).to(u.m))
    # interpolate total pressure, temperature and density profiles
    pdata_f = UnivariateSpline(hdata,np.array(np.log(data['p'])),k=1, s=s)
    Tdata_f = UnivariateSpline(hdata,np.array(np.log(data['T'])),k=1, s=s)
    rdata_f = UnivariateSpline(hdata,np.array(np.log(data['rho'])),k=1, s=s)
    #s=0.0 to ensure all points are strictly used for ionisation state
    muofT_f = UnivariateSpline(hdata,np.array(np.log(data['mu'])),k=1, s=0.0)

    outdata = Table()
    outdata['Z'] = Z
    outdata['p'] = np.exp(pdata_f(Z.to(u.m))) * data['p'].unit
    outdata['T'] = np.exp(Tdata_f(Z.to(u.m))) * data['T'].unit
    outdata['rho'] = np.exp(rdata_f(Z.to(u.m))) * data['rho'].unit
    outdata['mu'] = np.exp(muofT_f(Z.to(u.m))) * u.one

    return outdata


#----------------------------------------------------------------------------
# a simpler exponential atmosphere to test Spruit's analytical result
#----------------------------------------------------------------------------
def get_spruit_hs(
                   Z,
                   model_pars,
                   physical_constants,
                   option_pars
                 ):
    """ photospheric values of pressure and density are taken from VAL3c.
        Four options are available to select Alfven speed along the flux tube
        axis to be:
        constant, increase as the square root of Z, increase linearly and
        increase as the square 0f Z. We apply Bz~exp(-2z/chrom_scale) hence
        for Alfven speed \sqrt(B^2/rho) constant rho~exp(-4z/chrom_scale)...
        These are approximate due to the effect on density of the non-zero
        magnetic tension force.
        For HS equilibrium dp/dz = rho g., so cannot be isothermal?
    """
    p0 = model_pars['p0']
    r0 = 2.727e-07 * u.g/u.cm**3
    g0 = physical_constants['gravity']
    if option_pars['l_const']:
        pressure_Z = p0 * model_pars['chrom_scale']**3 /\
                            (model_pars['chrom_scale'] + Z)**3
        rho_Z = -p0 / g0 * 3. * model_pars['chrom_scale']**3/\
                            (model_pars['chrom_scale'] + Z)**4
        rtest = -p0 / g0 * 3. / model_pars['chrom_scale']
        model_pars['model'] = 'spruit_const'
    elif option_pars['l_sqrt']:
        pressure_Z = p0 *     model_pars['chrom_scale']**0.5/\
                             (model_pars['chrom_scale'] + Z)**0.5
        rho_Z = -0.5/g0 * p0 * model_pars['chrom_scale']**0.5/\
                             (model_pars['chrom_scale'] + Z)**1.5
        rtest = -0.5/g0 * p0 / model_pars['chrom_scale']
        model_pars['model'] = 'spruit_sqrt'
    elif option_pars['l_linear']:
        pressure_Z = p0 *     model_pars['chrom_scale']**1.5/\
                             (model_pars['chrom_scale'] + Z)**1.5
        rho_Z = -1.5/g0 * p0 * model_pars['chrom_scale']**1.5/\
                             (model_pars['chrom_scale'] + Z)**2.5
        rtest = -1.5/g0 * p0 / model_pars['chrom_scale']
        model_pars['model'] = 'spruit_linear'
    elif option_pars['l_square']:
        pressure_Z = p0 *     model_pars['chrom_scale']**3.5/\
                             (model_pars['chrom_scale'] + Z)**3.5
        rho_Z = -3.5/g0 * p0 * model_pars['chrom_scale']**3.5/\
                             (model_pars['chrom_scale'] + Z)**4.5
        rtest = -3.5/g0 * p0 / model_pars['chrom_scale']
        model_pars['model'] = 'spruit_square'
    else:
        raise ValueError("in hs_model.hs_atmosphere.get_spruit_hs set \
                  option_pars True for axial Alfven speed Z dependence")
    #to compare the derived density from hs-balance with VAL3c value:
    print('VAL rho(0) = ',r0.decompose(),' vs spruit rho(0) = ',rtest.decompose())
    Rgas_Z = u.Quantity(np.ones(Z.size), u.one)
    Rgas_Z *= physical_constants['boltzmann']/\
                physical_constants['proton_mass']/physical_constants['mu']

    return pressure_Z, rho_Z, Rgas_Z

#============================================================================
# Construct 3D hydrostatic profiles and include the magneto adjustments
#============================================================================
def vertical_profile(Z,
                     table,
                     magp0,
                     physical_constants, dz
                    ):
    """Return the vertical profiles for thermal pressure and density in 1D.
       Integrate in reverse from the corona to the photosphere to remove
       sensitivity to larger chromospheric gradients."""
    g0 = physical_constants['gravity'].to('m s-2')
    Rgas = u.Quantity(np.ones(table['Z'].size), u.one)
    Rgas *= (physical_constants['boltzmann']/\
                physical_constants['proton_mass']/table['mu']).to('m2 K-1 s-2')
    Rgas_Z  = Rgas[4:-4].copy()
    rdata   = u.Quantity(table['rho'], copy=True).to('kg m-3')
    rdata_Z = rdata[4:-4].copy()
    magp = magp0.to('kg m-1 s-2')
    # inverted SAC 4th order derivative scheme to minimise numerical error
    """evaluate upper boundary pressure from equation of state + enhancement,
       magp, which will be replaced by the mean magnetic pressure in the
       corona, then integrate from inner next pressure
    """
    table_T = u.Quantity(table['T'])
    linp_1 = table_T[-1]*rdata[-1]*Rgas[-1] + magp[-1]
    linp = u.Quantity(np.ones(len(Z)), unit=linp_1.unit)
    linp[-1] = table_T[-5]*rdata[-5]*Rgas[-5] + magp[-1]

#    for i in range(1,Z.size):
#        linp[-i-1] = (144.*linp[-i]+18.*linp[-i+1]
#                  -102.*(g0*rdata[-i-4]  )*dz
#                  - 84.*(g0*rdata[-i-5])*dz
#                  +  6.*(g0*rdata[-i-6])*dz
#                  )/162. + magp[-i-1] - magp[-i-0]
    for i in range(1,Z.size):
        linp[-i-1] = (1152.*linp[-i]
                  + 35.*(g0*rdata[-i-7])*dz
                  -112.*(g0*rdata[-i-6])*dz
                  -384.*(g0*rdata[-i-5])*dz
                  -784.*(g0*rdata[-i-4])*dz
                  + 77.*(g0*rdata[-i-3])*dz
                  )/1152. + magp[-i-1] - magp[-i-0]
    thermalp_Z = linp
    return thermalp_Z, rdata_Z, Rgas_Z


hmi_model = {'photo_scale': 0.6*u.Mm,       #scale height for photosphere
             'chrom_scale': 0.1*u.Mm,      #scale height for chromosphere
             'corona_scale': 2.5e3*u.Mm,      #scale height for the corona
             'coratio': 0.03*u.one,  #footpoint portion scaling as corona 
             'model': 'hmi_model',
             'phratio': 0.15*u.one,  #footpoint portion scaling as photosphere
             'pixel': 0.36562475*u.Mm,      #(HMI pixel)
             'radial_scale': 0.10979002*u.Mm, #=> FWHM = half pixel
             'nftubes': 1,
             'B_corona': 0.*u.T,
             'pBplus': 1e-3*u.T}
hmi_model['chratio'] = 1*u.one - hmi_model['coratio'] - hmi_model['phratio']

mfe_setup = {'photo_scale': 0.60*u.Mm,
             'chrom_scale': 0.4*u.Mm,
             'corona_scale': 0.25*u.Mm,  #scale height for the corona
             'coratio': 0.0*u.one,
             'model': 'mfe_setup',
             'phratio': 0.0*u.one,
             'pixel': 0.36562475*u.Mm,  #(HMI pixel)
             'radial_scale': 0.03938*u.Mm,
             #'radial_scale': 0.14*u.Mm,
             'nftubes': 1,
             #'B_corona': 4.85e-4*u.T,
             'B_corona': 5.5e-4*u.T,
             'pBplus': 12.0e-4*u.T}
mfe_setup['chratio'] = 1*u.one - mfe_setup['coratio'] - mfe_setup['phratio']
#if 1D or 2D set unused dimensions to 0, and unrequired xyz limits to 1.
mfe_setup['Nxyz'] = [128,128,128] # 3D grid
mfe_setup['Nxyz'] = [129,129,128] # 3D grid
mfe_setup['xyz']  = [-1*u.Mm,1*u.Mm,-1*u.Mm,1*u.Mm,3.6641221e-2*u.Mm,1.5877863*u.Mm] #grid size

spruit = {'photo_scale': 1.5*u.Mm,
          'chrom_scale': 0.5*u.Mm,
          'corona_scale': 100*u.Mm,      #scale height for the corona
          'coratio': 0.0*u.one,
          'model': 'spruit',
          'phratio': 0.0*u.one,
          'pixel': 0.1*u.Mm,              
          'radial_scale': 0.075*u.Mm,
          'nftubes': 1,
          'p0': 117200.0 * u.dyne/u.cm**2,
          'B_corona': 0.*u.T,
          'pBplus': 4.250e-4*u.T}
spruit['chratio'] = 1*u.one - spruit['coratio'] - spruit['phratio']
spruit['Nxyz'] = [128,128,256] # 3D grid
spruit['xyz']  = [-1.27*u.Mm,1.27*u.Mm,-1.27*u.Mm,1.27*u.Mm,0.0*u.Mm,25.5*u.Mm] #grid size


paper1 = {'photo_scale': 0.6*u.Mm,
          'chrom_scale': 0.1*u.Mm,
          'corona_scale': 2.5e3*u.Mm,         #scale height for the corona
          'coratio': 0.03*u.one,
          'model': 'paper1',
          'phratio': 0.0*u.one,
          'pixel': 0.36562475*u.Mm,              #(HMI pixel)
          'radial_scale': 0.10979002*u.Mm,
          'nftubes': 1,
          'B_corona': 9.2e-4*u.T,
          'pBplus': 1e-3*u.T}
paper1['chratio'] = 1*u.one - paper1['coratio'] - paper1['phratio']
paper1['Nxyz'] = [128,128,432] # 3D grid
paper1['xyz']  = [-1.27*u.Mm,1.27*u.Mm,-1.27*u.Mm,1.27*u.Mm,0.*u.Mm,8.62*u.Mm] #grid size
# uncomment to produce comaparable data to mfe_setup
#paper1['Nxyz'] = [127,128,128] # 3D grid
#paper1['xyz']  = [-1*u.Mm,1*u.Mm,-1*u.Mm,1*u.Mm,3.5e-2*u.Mm,1.6*u.Mm] #grid size

paper2a = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 2.5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2a',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 4,
           'B_corona': 8.4e-4*u.T,
           'pBplus': 1e-3*u.T}
paper2a['chratio'] = 1*u.one - paper2a['coratio'] - paper2a['phratio']
paper2a['Nxyz'] = [160,80,432] # 3D grid
paper2a['xyz']  = [-1.59*u.Mm,1.59*u.Mm,-0.79*u.Mm,0.79*u.Mm,0.*u.Mm,8.62*u.Mm] #grid size

paper2b = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 2.5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2b',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 4,
           'B_corona': 8.2e-4*u.T,
           'pBplus': 1.0e-3*u.T}
paper2b['chratio'] = 1*u.one - paper2b['coratio'] - paper2b['phratio']
paper2b['Nxyz'] = [50,50,140] # 3D grid
paper2b['xyz']  = [-0.49*u.Mm,0.49*u.Mm,-0.49*u.Mm,0.49*u.Mm,0*u.Mm,2.78*u.Mm] #grid size

paper2c = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2c',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 15,
           'B_corona': 5.95e-4*u.T,
           'pBplus': 1.0e-3*u.T}
paper2c['chratio'] = 1*u.one - paper2c['coratio'] - paper2c['phratio']
paper2c['Nxyz'] = [224,224,140] # 3D grid
paper2c['xyz']  = [-2.23*u.Mm,2.23*u.Mm,-2.23*u.Mm,2.23*u.Mm,0*u.Mm,2.78*u.Mm] #grid size

paper2d = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2d',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 18,
           'B_corona': 5.95e-4*u.T,
           'pBplus': 1.0e-3*u.T}
paper2d['chratio'] = 1*u.one - paper2d['coratio'] - paper2d['phratio']
paper2d['Nxyz'] = [224,224,140] # 3D grid
paper2d['xyz']  = [-2.23*u.Mm,2.23*u.Mm,-0.79*u.Mm,0.79*u.Mm,0*u.Mm,2.78*u.Mm] #grid size


def get_coords(Nxyz, xyz):
    """
    get_coords returns a non-dimensional dictionary describing the domain
    coordinates.
    """
    dz=(xyz[5]-xyz[4])/(Nxyz[2]-1)
    Z    = u.Quantity(np.linspace(xyz[4].value, xyz[5].value, Nxyz[2]), unit=xyz.unit)
    Zext = u.Quantity(np.linspace(Z.min().value-4.*dz.value, Z.max().value+4.*dz.value, Nxyz[2]+8), unit=Z.unit)
    coords = {
              'dx':(xyz[1]-xyz[0])/(Nxyz[0]-1),
              'dy':(xyz[3]-xyz[2])/(Nxyz[1]-1),
              'dz':(xyz[5]-xyz[4])/(Nxyz[2]-1),
              'xmin':xyz[0],
              'xmax':xyz[1],
              'ymin':xyz[2],
              'ymax':xyz[3],
              'zmin':xyz[4],
              'zmax':xyz[5],
              'Z':Z,
              'Zext':Zext
             }

    return coords
#-----------------------------------------------------------------------------
#
def get_hmi_map(
                model_pars, option_pars,
                indx, 
                dataset = 'hmi_m_45s_2014_07_06_00_00_45_tai_magnetogram_fits', 
                sunpydir = os.path.expanduser('~/sunpy/data/'),
                figsdir = os.path.expanduser('~/figs/hmi/'),
                l_newdata = False
               ):
    """ indx is 4 integers: lower and upper indices each of x,y coordinates 
#    dataset of the form 'hmi_m_45s_2014_07_06_00_00_45_tai_magnetogram_fits'
#    """
    from scipy.interpolate import RectBivariateSpline
    from sunpy.net import vso
    import sunpy.map
    client = vso.VSOClient()
    results = client.query(vso.attrs.Time("2014/07/05 23:59:50",
                                          "2014/07/05 23:59:55"), 
                           vso.attrs.Instrument('HMI'),
                           vso.attrs.Physobs('LOS_magnetic_field'))
    #print results.show()                       

    if l_newdata:
        if not os.path.exists(sunpydir):
            raise ValueError("in get_hmi_map set 'sunpy' dir for vso data\n"+ 
        "for large files you may want link to local drive rather than network")
        client.get(results).wait(progress=True)
    if not os.path.exists(figsdir):
        os.makedirs(figsdir)

    hmi_map = sunpy.map.Map(sunpydir+dataset)
    #hmi_map = hmi_map.rotate()
    #hmi_map.peek()
    s = hmi_map.data[indx[0]:indx[1],indx[2]:indx[3]] #units of Gauss Bz
    s *= u.G
    nx = s.shape[0]
    ny = s.shape[1]
    nx2, ny2 = 2*nx, 2*ny # size of interpolant 
    #pixel size in arc seconds
    dx, dy = hmi_map.scale.items()[0][1],hmi_map.scale.items()[1][1] 
    x, y = np.mgrid[
             hmi_map.xrange[0]+indx[0]*dx:hmi_map.xrange[0]+indx[1]*dx:1j*nx2,
             hmi_map.xrange[0]+indx[2]*dy:hmi_map.xrange[0]+indx[3]*dy:1j*ny2
                     ]
    #arrays to interpolate s from/to
    fx =   u.Quantity(np.linspace(x.min().value,x.max().value,nx), unit=x.unit)
    fy =   u.Quantity(np.linspace(y.min().value,y.max().value,ny), unit=y.unit)
    xnew = u.Quantity(np.linspace(x.min().value,x.max().value,nx2), unit=x.unit)
    ynew = u.Quantity(np.linspace(y.min().value,y.max().value,ny2), unit=y.unit)
    f  = RectBivariateSpline(fx,fy,s.to(u.T))
    #The initial model assumes a relatively small region, so a linear 
    #Cartesian map is applied here. Consideration may be required if larger
    #regions are of interest, where curvature or orientation near the lim
    #of the surface is significant. 
    s_int  = f(xnew,ynew) #interpolate s and convert units to Tesla
    s_int /= 4. # rescale s as extra pixels will sum over FWHM
    x_int  = x  * 7.25e5 * u.m    #convert units to metres
    y_int  = y  * 7.25e5 * u.m
    dx_int = dx * 7.25e5 * u.m
    dy_int = dy * 7.25e5 * u.m 
    FWHM  = 0.5*(dx_SI+dy_SI)
    smax  = max(abs(s.min()),abs(s.max())) # set symmetric plot scale
    cmin  = -smax*1e-4
    cmax  =  smax*1e-4
#    
#    filename = 'hmi_map'
#    import loop_plots as mhs
#    mhs.plot_hmi(
#             s*1e-4,x_SI.min(),x_SI.max(),y_SI.min(),y_SI.max(),
#             cmin,cmax,filename,savedir,annotate = '(a)'
#            )
#    filename = 'hmi_2x2_map'
#    mhs.plot_hmi(
#             s_SI*4,x_SI.min(),x_SI.max(),y_SI.min(),y_SI.max(),
#             cmin,cmax,filename,savedir,annotate = '(a)'
#            )
#
#    return s_SI, x_SI, y_SI, nx2, ny2, dx_SI, dy_SI, cmin, cmax, FWHM
    dz=(xyz[5]-xyz[4])/(Nxyz[2]-1)
    Z    = u.Quantity(np.linspace(xyz[4].value, xyz[5].value, Nxyz[2]), unit=xyz.unit)
    Zext = u.Quantity(np.linspace(Z.min().value-4.*dz.value, Z.max().value+4.*dz.value, Nxyz[2]+8), unit=Z.unit)
    coords = {
              'dx':(xyz[1]-xyz[0])/(Nxyz[0]-1),
              'dy':(xyz[3]-xyz[2])/(Nxyz[1]-1),
              'dz':(xyz[5]-xyz[4])/(Nxyz[2]-1),
              'xmin':xyz[0],
              'xmax':xyz[1],
              'ymin':xyz[2],
              'ymax':xyz[3],
              'zmin':xyz[4],
              'zmax':xyz[5],
              'Z':Z,
              'Zext':Zext
             }

    return coords

##============================================================================
##option parameters
##============================================================================
def set_options(model, l_mpi, l_gdf=True):

    """This module assigns the logical options for the model. If adding 
    new models with additional logical arguments add it to the default 
    list as false, include an if statement for True update the
    dictionary option_pars 
    
    """    
    #default arguments
    option_pars = {
        'l_hdonly':      False,# set mag field zero to check background
        'l_ambB':        False,# include some ambient magnetic field b_z
        'l_spruit':      False,# thin flux tube model to check Spruit 
        'l_const':       False,# axial Alfven speed const  Z-depend (Spruit)
        'l_sqrt':        False,# axial Alfven speed sqrt   Z-depend (Spruit)
        'l_linear':      False,# axial Alfven speed linear Z-depend (Spruit)
        'l_square':      False,# axial Alfven speed square Z-depend (Spruit)
        'l_B0_expz':     False,# Z-depend of Bz(r=0) exponentials
        'l_B0_quadz':    False,# Z-depend of Bz(r=0) polynomials + exponential 
        'l_B0_rootz':    False,# Z-depend of Bz(r=0) sqrt polynomials 
        'l_single':      False,# only one flux tube
        'l_hmi':         False,# construct photopheric map of Bz from HMI/SDI
        'l_tube_pair':   False,# pair of flux tubes
        'l_multi_netwk': False,# multiple flux tubes as described in GFE (2014)
        'l_multi_lanes': False,# multiple flux tubes as described in GFE (2014)
        'l_multi_twist': False,# multiple flux tubes as described in GFE (2014)
        'l_2D_loop':     False,# make a 2D loop with sinusoidal Bz(x,0,0)
        'l_mfe':         False,# model Viktor's model from MFE (2014)      
        'l_atmos_val3c_mtw':False,# interpolate composite VAL3c+MTW atmosphere
        'suffix':        '.gdf'
    }
    #revise optional parameters depending on configuration required
    if model['model'] == 'hmi_model':
        option_pars['l_hmi']             = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'mfe_setup':
        option_pars['l_single']          = True 
        option_pars['l_mfe']             = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'spruit':    
        option_pars['l_single']          = True 
        option_pars['l_spruit']          = True 
    if model['model'] == 'paper1':
        option_pars['l_ambB']            = True 
        option_pars['l_B0_expz']         = True
        option_pars['l_single']          = True
        option_pars['l_atmos_val3c_mtw'] = True
    if model['model'] == 'paper2a':
        option_pars['l_ambB']            = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_tube_pair']       = True 
        option_pars['l_atmos_val3c_mtw'] = True
    if model['model'] == 'paper2b':
        option_pars['l_ambB']            = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_multi_twist'    ] = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'paper2c':
        option_pars['l_ambB']            = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_multi_netwk']     = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'paper2d':
        option_pars['l_ambB']            = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_multi_lanes'    ] = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'hmi_model':
        option_pars['l_B0_quadz']        = True 
        option_pars['l_single']          = True 
        option_pars['l_hmi']             = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'loop_model':
        option_pars['l_B0_quadz']        = True 
        option_pars['l_single']          = True 
        option_pars['l_2D_loop']         = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if l_mpi:
        option_pars['l_mpi'] = True
    else:
        option_pars['l_mpi'] = False
    return option_pars

def get_parameters():
#============================================================================
    # Dimensional units in terms of SI
#============================================================================
    scales   = {
            'length':         1e2*u.Mm,
            'density':        1e-4*u.kg/u.m**3,
            'velocity':       1e2*u.m/u.s,
            'temperature':    1.0*u.K, 
            'magnetic':       1e-3*u.T #mT
           }
    scales['energy density'] = scales['density'] * scales['velocity']**2
    scales['time'] = scales['length'] / scales['velocity'] 
    scales['mass'] = scales['density'] * scales['length']**3 
    scales['force density'] = scales['density'] * scales['velocity'] / \
                       scales['time'] #D momentum/Dt force density balance 
#============================================================================
# physical constants
#============================================================================
    physical_constants = {'gamma':      5.0/3.0         , 
                          'mu':         0.602           , 
                          'mu0':        asc.mu0         , 
                          'boltzmann':  asc.k_B         ,
                          'proton_mass':asc.m_p         ,
                          'gravity':    -2.74e2*u.m/u.s/u.s
                         }

    return scales, physical_constants
