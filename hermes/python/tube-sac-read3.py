import numpy as np
import vtk
from scipy import special
import scipy.io
import scipy.ndimage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import use
#use('Agg') 
from matplotlib import pyplot
import sip

from pylab import *
import struct

from mayavi import mlab

sip.setapi('QString', 2)
from vtk.util import numpy_support as VN

max_slice=0


file = open('../../../temp/washmc__480000.out','rb')



file.seek(0,2)
eof = file.tell()
file.seek(0,0)

name = file.read(79)

nit = fromfile(file,dtype=int32,count=1)

t = fromfile(file,dtype=float64,count=1)
ndim=fromfile(file,dtype=int32,count=1)
neqpar=fromfile(file,dtype=int32,count=1)
nw=fromfile(file,dtype=int32,count=1)

ndata = fromfile(file,dtype=int32,count=ndim)[:ndim]

varbuf = fromfile(file,dtype=float,count=7)[:7]

#if ndim=2
varnames = file.read(79)

#if ndim=3
 
 
 
 
 
 #typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b} CEV;

if ndim==3:
    alldat=fromfile(file,dtype=float,count=(nw+ndim)*ndata[0]*ndata[1]*ndata[2])[:(nw+ndim)*ndata[0]*ndata[1]*ndata[2]]
    #if size(alldat)<(nw+ndim)*ndata[0]*ndata[1]*ndata[2]:
    #    alldat=resize(alldat,(nw+ndim)*ndata[0]*ndata[1]*ndata[2])
    alldat=np.reshape(alldat,(nw+ndim,ndata[0],ndata[1],ndata[2],),'C')

file.close()





#x = np.zeros(ndata[0]*ndata[1]*ndata[2])
#y = np.zeros(ndata[0]*ndata[1]*ndata[2])
#z = np.zeros(ndata[0]*ndata[1]*ndata[2])

#for i in range(data.GetNumberOfPoints()):
#        x[i],y[i],z[i] = data.GetPoint(i)
        
        

#x = x.reshape(ndata,order='F')
#y = y.reshape(ndata,order='F')
#z = z.reshape(ndata,order='F')

x=alldat[0,:,:,:]
y=alldat[1,:,:,:]
z=alldat[2,:,:,:]

#Bx=alldat[13,:,:,:]+alldat[5,:,:,:]
#By=alldat[14,:,:,:]+alldat[6,:,:,:]
#Bz=alldat[15,:,:,:]+alldat[7,:,:,:]

Bx=alldat[5,:,:,:]
By=alldat[6,:,:,:]
Bz=alldat[7,:,:,:]


dens=alldat[0,:,:,:]

scalar = mlab.pipeline.scalar_field(Bz)
vector = mlab.pipeline.vector_field(Bx, By, Bz)

densscalar = mlab.pipeline.scalar_field(dens)

#### Visualize the field ####################################################

fig = mlab.figure(1, size=(800, 800), bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))

# Vector field first for seeing the field lines

magnitude = mlab.pipeline.extract_vector_norm(vector)
field_lines = mlab.pipeline.streamline(magnitude)

# Scalar field for visualisation image planes 
# First the photosphere
mlab.pipeline.image_plane_widget(densscalar, plane_orientation='z_axes', slice_index=0, colormap='RdBu', transparent=False, opacity=1)

# Upper layers

for i in range(0, max_slice):
    
    mlab.pipeline.image_plane_widget(densscalar, plane_orientation='z_axes', slice_index=50, colormap='RdBu', transparent=True, opacity=0.5)




engine=mlab.get_engine()
image_plane_widget = engine.scenes[0].children[2].children[0].children[0]
image_plane_widget.ipw.origin = array([ 0.5,  0.5,  4. ])
image_plane_widget.ipw.origin = array([ 0.5,  0.5,  4. ])
image_plane_widget.ipw.slice_index = 3
image_plane_widget.ipw.point1 = array([ 128.5,    0.5,    4. ])
image_plane_widget.ipw.point1 = array([ 128.5,    0.5,    4. ])
image_plane_widget.ipw.point2 = array([   0.5,  128.5,    4. ])
image_plane_widget.ipw.point2 = array([   0.5,  128.5,    4. ])
image_plane_widget.ipw.slice_position = 4.0

streamlines = engine.scenes[0].children[1].children[0].children[0].children[0]


streamlines.seed.widget = streamlines.seed.widget_list[2]
streamlines.seed.widget.enabled = False
streamlines.seed.widget.interactor = None
# streamlines.seed.widget = <tvtk.tvtk_classes.plane_widget.PlaneWidget object at 0x136851d70>

streamlines.seed.widget.enabled = False

streamlines.seed.widget.point1 = array([ 64.5 ,  96.25,  32.75])
streamlines.seed.widget.origin = array([ 32.75,  32.75,  64.5 ])
streamlines.seed.widget.center = array([ 64.5,  64.5,  64.5])
streamlines.seed.widget.normal = array([ 0.,  0.,  1.])

streamlines.seed.widget.point2 = array([ 32.75,  96.25,  64.5 ])
streamlines.seed.widget.normal_to_z_axis = True
streamlines.stream_tracer.start_position = array([ 0.,  0.,  0.])
streamlines.stream_tracer.progress = 1.0
streamlines.stream_tracer.integration_direction = 'both'
streamlines.stream_tracer.start_position = array([ 0.,  0.,  0.])
streamlines.stream_tracer.maximum_propagation = 500.0
# streamlines.clean_filter.input_connection = <tvtk.tvtk_classes.algorithm_output.AlgorithmOutput object at 0x143d3c890>
streamlines.tube_filter.default_normal = array([ 0.,  0.,  1.])
# streamlines.tube_filter.input_connection = <tvtk.tvtk_classes.algorithm_output.AlgorithmOutput object at 0x143d3c890>
streamlines.actor.mapper.progress = 1.0
streamlines.actor.mapper.scalar_range = array([ 0.        ,  0.09823792])
# streamlines.actor.mapper.input = <tvtk.tvtk_classes.poly_data.PolyData object at 0x143d3c890>
streamlines.streamline_type = 'tube'

streamlines.tube_filter.progress = 1.0
streamlines.tube_filter.vary_radius = 'vary_radius_by_scalar'
streamlines.tube_filter.number_of_sides = 6
streamlines.tube_filter.default_normal = array([ 0.,  0.,  1.])
streamlines.tube_filter.capping = True



streamlines.actor.mapper.scalar_range = array([ 0.        ,  0.09823792])
# streamline.actor.mapper.input = <tvtk.tvtk_classes.poly_data.PolyData object at 0x1416a2230>
streamlines.streamline_type = 'line'

streamlines.seed.widget = streamlines.seed.widget_list[0]
streamlines.seed.widget.radius = 10.0
streamlines.seed.widget.center = array([  64.43284748,   64.38386692,  110.6832856 ])
streamlines.seed.widget.handle_direction = array([ 1.,  0.,  0.])

# streamline.seed.widget = <tvtk.tvtk_classes.sphere_widget.SphereWidget object at 0x138199890>
streamlines.seed.widget.enabled = False

streamlines.seed.widget.enabled = True
streamlines.seed.widget.handle_direction = array([ 1.,  0.,  0.])

streamlines.seed.widget.handle_direction = array([ 1.,  0.,  0.])
streamlines.seed.widget.handle_direction = array([ 1.,  0.,  0.])


streamlines.seed.widget.enabled = False


#streamlines.



#mlab.colorbar(magnitude)
#mlab.title('polar mesh')
#mlab.outline(magnitude)
mlab.axes(magnitude)


scene = engine.scenes[0]
scene.scene.save(u'/Users/mike/proj/solar/hermes/python/snapshot.png')
#field_lines = mlab.pipeline.streamline(magnitude, seedtype='line', integration_direction='both', colormap='bone', vmin=0, vmax=1)

#streamlines.stream_tracer.maximum_propagation = 100.
#streamlines.seed.widget.point1 = [69, 75.5, 75.5]
#streamlines.seed.widget.point2 = [82, 75.5, 75.5]
#streamlines.seed.widget.resolution = 50
#streamlines.seed.widget.enabled = False

#mlab.view(42, 73, 104, [79,  75,  76])

mlab.show_pipeline()
mlab.show()
