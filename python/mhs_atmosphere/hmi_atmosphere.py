# -*- coding: utf-8 -*-
"""
Created on Wed 27th May 14:36:00 2015

@author: sm1fg

This is the main module to construct a magnetohydrostatic solar atmosphere,
given a magnetic network of self-similar magnetic flux tubes based on an HMI
data set of the line of sight magnetic field component and 
save the output to gdf format.

To select an existing configuration change the import as model_pars, set Nxyz,
xyz_SI and any other special parameters, then execute mhs_atmopshere.

To add new configurations:
add the model options to set_options in parameters/options.py;
add options required in parameters/model_pars.py;
add alternative empirical data sets to hs_model/;
add alternativ table than interploate_atmosphere in hs_model/hs_atmosphere.py;
add option to get_flux_tubes in mhs_model/flux_tubes.py

If an alternative formulation of the flux tube is required add options to
construct_magnetic_field and construct_pairwise_field in
mhs_model/flux_tubes.py

Plotting options are included in plot/mhs_plot.py
"""

import os
import numpy as np
import pysac.mhs_atmosphere as atm
import astropy.units as u
from pysac.mhs_atmosphere.parameters.model_pars import hmi_model as model_pars
import astropy
l_astropy_0=True
if astropy.__version__[0]=='1':
    l_astropy_0=False
#==============================================================================
#check whether mpi is required and the number of procs = size
#==============================================================================
try:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    collective=True
    l_mpi = True
    l_mpi = l_mpi and (size != 1)
except ImportError:
    collective=False
    l_mpi = False
    rank = 0
    size = 1
#==============================================================================
#set up model parameters
#==============================================================================
local_procs=1
#standard set of logical switches
option_pars = atm.options.set_options(model_pars, l_mpi, l_gdf=True)
#standard conversion to dimensionless units and physical constants
scales, physical_constants = \
    atm.units_const.get_parameters()

if not l_astropy_0:
    dataset = 'hmi_m_45s_2014_07_06_00_00_45_tai_magnetogram_fits'
else:
    dataset = 'hmi_m_45s_2014_07_06_00_00_45_tai_magnetogram.fits'
l_newdata = True # change this to False if a local copy already exists at ~/sunpy/data/
#obtain code coordinates and model parameters in astropy units
model_pars['Nxyz'] = [64,64,64] # 3D grid
model_pars['xyz']  = [-0.63*u.Mm,0.63*u.Mm,-0.63*u.Mm,0.63*u.Mm,0.0*u.Mm,3.8*u.Mm]
indx = [1785,1800,1814,1829]
interpfactor=5
frac = [0.25,0.25]
coords = atm.model_pars.get_hmi_coords(
                model_pars['Nxyz'],
                model_pars['xyz'],
                indx = indx,
                dataset = dataset,
                l_newdata = False,
                interpfactor=interpfactor,
                frac=frac,
                rank=rank,
                lmpi=l_mpi
               )
model_pars['xyz'][0:4] = coords['xmin'],coords['xmax'],coords['ymin'],coords['ymax']

#interpolate the hs 1D profiles from empirical data source[s]
empirical_data = atm.hs_atmosphere.read_VAL3c_MTW(mu=physical_constants['mu'])

table = \
        atm.hs_atmosphere.interpolate_atmosphere(empirical_data,
                                   coords['Zext']
                                  )

#==============================================================================
#calculate 1d hydrostatic balance from empirical density profile
#==============================================================================
# the hs pressure balance is enhanced by pressure equivalent to the
# residual mean coronal magnetic pressure contribution once the magnetic
# field has been applied
magp_meanz = np.ones(len(coords['Z'])) * u.one
magp_meanz *= model_pars['pBplus']**2/(2*physical_constants['mu0'])

pressure_Z, rho_Z, Rgas_Z = atm.hs_atmosphere.vertical_profile(
                                                 coords['Z'],
                                                 table,
                                                 magp_meanz,
                                                 physical_constants,
                                                 coords['dz']
                                                 )

#==============================================================================
# load flux tube footpoint parameters
#==============================================================================
# axial location and value of Bz at each footpoint
Stmp,xtmp,ytmp,FWHM,sdummy,xdummy,ydummy,cmax,cmin \
               =atm.parameters.model_pars.get_hmi_map(
                indx,
                dataset = dataset,
                l_newdata = False,
                interpfactor=interpfactor,
                frac=frac,
                rank=rank,
                lmpi=l_mpi
               )
xi, yi, Si = xtmp.to(u.Mm).reshape(xtmp.size),\
             ytmp.to(u.Mm).reshape(ytmp.size), Stmp.reshape(Stmp.size)
model_pars['radial_scale'] = 0.5*FWHM.to(u.Mm)/np.sqrt(2*np.log(2))
#==============================================================================
# split domain into processes if mpi
#==============================================================================
ax, ay, az = np.mgrid[coords['xmin']:coords['xmax']:1j*model_pars['Nxyz'][0],
                      coords['ymin']:coords['ymax']:1j*model_pars['Nxyz'][1],
                      coords['zmin']:coords['zmax']:1j*model_pars['Nxyz'][2]]

axindex=np.arange(model_pars['Nxyz'][0])
# split the grid between processes for mpi
if l_mpi:
    x_chunks = np.array_split(ax, size, axis=0)
    y_chunks = np.array_split(ay, size, axis=0)
    z_chunks = np.array_split(az, size, axis=0)
    i_chunks = np.array_split(axindex, size, axis=0)

    x = comm.scatter(x_chunks, root=0)
    y = comm.scatter(y_chunks, root=0)
    z = comm.scatter(z_chunks, root=0)
    xindex=i_chunks[rank]
else:
    x, y, z = ax, ay, az
    xindex=axindex

x = u.Quantity(x, unit=coords['xmin'].unit)
y = u.Quantity(y, unit=coords['ymin'].unit)
z = u.Quantity(z, unit=coords['zmin'].unit)
#==============================================================================
# initialize zero arrays in which to add magnetic field and mhs adjustments
#==============================================================================
Bx   = u.Quantity(np.zeros(x.shape), unit=u.T)  # magnetic x-component
By   = u.Quantity(np.zeros(x.shape), unit=u.T)  # magnetic y-component
Bz   = u.Quantity(np.zeros(x.shape), unit=u.T)  # magnetic z-component
pressure_m = u.Quantity(np.zeros(x.shape), unit=u.Pa) # magneto-hydrostatic adjustment to pressure
rho_m = u.Quantity(np.zeros(x.shape), unit=u.kg/u.m**3)      # magneto-hydrostatic adjustment to density
# initialize zero arrays in which to add balancing forces and magnetic tension
Fx   = u.Quantity(np.zeros(x.shape), unit=u.N/u.m**3)  # balancing force x-component
Fy   = u.Quantity(np.zeros(x.shape), unit=u.N/u.m**3)  # balancing force y-component
# total tension force for comparison with residual balancing force
Btensx  = u.Quantity(np.zeros(x.shape), unit=u.N/u.m**3)
Btensy  = u.Quantity(np.zeros(x.shape), unit=u.N/u.m**3)
#==============================================================================
#calculate the magnetic field and pressure/density balancing expressions
#==============================================================================
for i in range(0,Si.size):
    for j in range(i,Si.size):
        if rank == 0:
            print'calculating ij-pair:',i,j
        if i == j:
            pressure_mi, rho_mi, Bxi, Byi ,Bzi, B2x, B2y =\
                atm.flux_tubes.construct_magnetic_field(
                                             x, y, z,
                                             xi[i], yi[i], Si[i],
                                             model_pars, option_pars,
                                             physical_constants,
                                             scales
                                            )
            Bx, By, Bz = Bxi+Bx, Byi+By ,Bzi+Bz
            Btensx += B2x
            Btensy += B2y
            pressure_m += pressure_mi
            rho_m += rho_mi
        else:
            pressure_mi, rho_mi, Fxi, Fyi, B2x, B2y =\
                atm.flux_tubes.construct_pairwise_field(
                                             x, y, z,
                                             xi[i], yi[i],
                                             xi[j], yi[j], Si[i], Si[j],
                                             model_pars,
                                             option_pars,
                                             physical_constants,
                                             scales
                                            )
            pressure_m += pressure_mi
            rho_m += rho_mi
            Fx   += Fxi
            Fy   += Fyi
            Btensx += B2x
            Btensy += B2y

# clear some memory
del pressure_mi, rho_mi, Bxi, Byi ,Bzi, B2x, B2y

import matplotlib.pyplot as plt
plt.figure()
plt.pcolormesh(x[:,:,0].T.value,y[:,:,0].T.value,Bz[:,:,0].T.value,vmin=cmin,vmax=cmax)
plt.xlabel('lon [Mm]')
plt.ylabel('lat [Mm]')
plt.axis([x.min().value,x.max().value,y.min().value,y.max().value])
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'$B_z$ [T]')
cbar.solids.set_edgecolor("face")
plt.savefig('model_derived_hmi.png')
plt.close()
print 'Bz[:,:,0].min(), Bz[:,:,0].max()=', Bz[:,:,0].min(), Bz[:,:,0].max()
#==============================================================================
# Construct 3D hs arrays and then add the mhs adjustments to obtain atmosphere
#==============================================================================
# select the 1D array spanning the local mpi process; the add/sub of dz to
# ensure all indices are used, but only once
indz = np.where(coords['Z'] >= z.min()-0.1*coords['dz']) and \
       np.where(coords['Z'] <= z.max()+0.1*coords['dz'])
pressure_z, rho_z, Rgas_z = pressure_Z[indz], rho_Z[indz], Rgas_Z[indz]
# local proc 3D mhs arrays
pressure, rho = atm.mhs_3D.mhs_3D_profile(z,
                                   pressure_z,
                                   rho_z,
                                   pressure_m,
                                   rho_m
                                  )
magp = (Bx**2 + By**2 + Bz**2)/(2.*physical_constants['mu0'])
if rank ==0:
    print'max pB corona = ',magp[:,:,-1].max().decompose()
    print'min pB corona = ',magp[:,:,-1].min().decompose()
energy = atm.mhs_3D.get_internal_energy(pressure,
                                                  magp,
                                                  physical_constants)
Rgas = u.Quantity(np.zeros(x.shape), unit=Rgas_z.unit)
Rgas[:] = Rgas_z
temperature = pressure/rho/Rgas
if not option_pars['l_hdonly']:
    inan = np.where(magp <=1e-7*pressure.min())
    magpbeta = magp
    magpbeta[inan] = 1e-7*pressure.min()  # low pressure floor to avoid NaN
    pbeta  = pressure/magpbeta
else:
    pbeta  = magp+1.0    #dummy to avoid NaN
alfven = np.sqrt(2.*magp/rho)
#if rank == 0:
#    print'Alfven speed Z.min to Z.max =',\
#    alfven[model_pars['Nxyz'][0]/2,model_pars['Nxyz'][1]/2, 0].decompose(),\
#    alfven[model_pars['Nxyz'][0]/2,model_pars['Nxyz'][1]/2,-1].decompose()
cspeed = np.sqrt(physical_constants['gamma']*pressure/rho)
#============================================================================
# Save data for SAC and plotting
#============================================================================
# set up data directory and file names
# may be worthwhile locating on /data if files are large
datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+model_pars['model']+'/')
filename = datadir + model_pars['model'] + option_pars['suffix']
if not os.path.exists(datadir):
    os.makedirs(datadir)
sourcefile = datadir + model_pars['model'] + '_sources' + option_pars['suffix']
aux3D = datadir + model_pars['model'] + '_3Daux' + option_pars['suffix']
aux1D = datadir + model_pars['model'] + '_1Daux' + option_pars['suffix']

# save the variables for the initialisation of a SAC simulation
atm.mhs_snapshot.save_SACvariables(
          filename,
          rho,
          Bx,
          By,
          Bz,
          energy,
          option_pars,
          physical_constants,
          coords,
          model_pars['Nxyz'],
          xindex,
          rank=rank,
          collective=collective
         )
# save the balancing forces as the background source terms for SAC simulation
atm.mhs_snapshot.save_SACsources(
          sourcefile,
          Fx,
          Fy,
          option_pars,
          physical_constants,
          coords,
          model_pars['Nxyz'],
          xindex,
          rank=rank,
          collective=collective
         )
# save auxilliary variable and 1D profiles for plotting and analysis
atm.mhs_snapshot.save_auxilliary3D(
          aux3D,
          pressure_m,
          rho_m,
          temperature,
          pbeta,
          alfven,
          cspeed,
          Btensx,
          Btensy,
          option_pars,
          physical_constants,
          coords,
          model_pars['Nxyz'],
          xindex,
          rank=rank,
          collective=collective
         )
if rank==0:
    atm.mhs_snapshot.save_auxilliary1D(
          aux1D,
          pressure_Z,
          rho_Z,
          Rgas_Z,
          option_pars,
          physical_constants,
          coords,
          model_pars['Nxyz'],
          rank=rank,
          collective=False
         )
if rho.min()<0 or pressure.min()<0:
    print"FAIL rank {}: negative rho.min() {} and/or pressure.min() {}.".format(
                            rank,rho.min(),pressure.min())
FWHM = 2*np.sqrt(np.log(2))*model_pars['radial_scale']
print'FWHM(0) =',FWHM
