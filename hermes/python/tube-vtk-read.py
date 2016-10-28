import numpy as np
import vtk
from scipy import special
import scipy.io
import scipy.ndimage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sip

from mayavi import mlab

sip.setapi('QString', 2)
from vtk.util import numpy_support as VN

max_slice=1

#reader = vtkStructuredGridReader()
reader = vtk.vtkDataSetReader()
reader.SetFileName("bb480.vtk")
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()

data = reader.GetOutput()

dim = data.GetDimensions()
vec = list(dim)
vec = [i-1 for i in dim]
vec.append(3)


#test=data.GetCellData().GetArray('bb')

#bb = VN.vtk_to_numpy(data.GetCellData().GetArray('bb'))
#bb = VN.vtk_to_numpy(data.GetStructuredGridOutput().GetData('bb'))

bbd=data.GetPointData().GetArray('bb')
bb=VN.vtk_to_numpy(bbd)
#b = u.reshape(vec,order='F')
bb = bb.reshape([128,128,128,3],order='F')

x = np.zeros(data.GetNumberOfPoints())
y = np.zeros(data.GetNumberOfPoints())
z = np.zeros(data.GetNumberOfPoints())

for i in range(data.GetNumberOfPoints()):
        x[i],y[i],z[i] = data.GetPoint(i)

x = x.reshape(dim,order='F')
y = y.reshape(dim,order='F')
z = z.reshape(dim,order='F')

Bx=bb[:,:,:,1]
By=bb[:,:,:,2]
Bz=bb[:,:,:,0]

scalar = mlab.pipeline.scalar_field(Bz)
vector = mlab.pipeline.vector_field(Bx, By, Bz)

#### Visualize the field ####################################################

fig = mlab.figure(1, size=(800, 800), bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))

# Vector field first for seeing the field lines

magnitude = mlab.pipeline.extract_vector_norm(vector)
field_lines = mlab.pipeline.streamline(magnitude)

# Scalar field for visualisation image planes 
# First the photosphere
mlab.pipeline.image_plane_widget(scalar, plane_orientation='z_axes', slice_index=0, colormap='RdBu', transparent=False, opacity=1)

# Upper layers

for i in range(0, max_slice):
    
    mlab.pipeline.image_plane_widget(scalar, plane_orientation='z_axes', slice_index=50, colormap='RdBu', transparent=True, opacity=0.5)



#mlab.colorbar(magnitude)
#mlab.title('polar mesh')
#mlab.outline(magnitude)
mlab.axes(magnitude)


#field_lines = mlab.pipeline.streamline(magnitude, seedtype='line', integration_direction='both', colormap='bone', vmin=0, vmax=1)

#field_lines.stream_tracer.maximum_propagation = 100.
#field_lines.seed.widget.point1 = [69, 75.5, 75.5]
#field_lines.seed.widget.point2 = [82, 75.5, 75.5]
#field_lines.seed.widget.resolution = 50
#field_lines.seed.widget.enabled = False

#mlab.view(42, 73, 104, [79,  75,  76])
mlab.show_pipeline()
mlab.show()