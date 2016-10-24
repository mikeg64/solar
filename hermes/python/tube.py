
#start using python tube.py
#before starting set the following

#export QT_API=pyqt
#export ETS_TOOLKIT=qt4

#(avoids error ValueError: API 'QString' has already been set to version 1 )

import numpy as np
from scipy import special
import scipy.io
import scipy.ndimage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sip

from mayavi import mlab

sip.setapi('QString', 2)

max_slice = 3

#### Read IDL save file ###################################################

idl_save = scipy.io.readsav('../../../temp/test.sav',python_dict=True, verbose=False)

#### Unpack variables #####################################################

date = idl_save['b3dbox']['date'][0]
jd = idl_save['b3dbox']['jd'][0]
x = idl_save['b3dbox']['x'][0]
y = idl_save['b3dbox']['y'][0]
z = idl_save['b3dbox']['z'][0]
l = idl_save['b3dbox']['l'][0]
b = idl_save['b3dbox']['b'][0]


Bx = idl_save['b3dbox']['bx'][0].T
By = idl_save['b3dbox']['by'][0].T
Bz = idl_save['b3dbox']['bz'][0].T

#### Vector field forced into box shape ###################################
# The new box size will be: x, x, x from x, y, z

Bx = scipy.ndimage.zoom(Bx,[1,1,3], order=0)
By = scipy.ndimage.zoom(By,[1,1,3], order=0)
Bz = scipy.ndimage.zoom(Bz,[1,1,3], order=0)

#### Vector field and scalar field creation ################################
# Scalar field = only Bz component!!!

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

