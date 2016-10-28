import numpy
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN

reader = vtkStructuredPointsReader()
reader.SetFileName('../../../temp/bb480.vtk')
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()

data = reader.GetOutput()

dim = data.GetDimensions()
vec = list(dim)
vec = [i-1 for i in dim]
vec.append(3)

u = VN.vtk_to_numpy(data.GetCellData().GetArray('velocity'))
b = VN.vtk_to_numpy(data.GetCellData().GetArray('cell_centered_B'))

u = u.reshape(vec,order='F')
b = b.reshape(vec,order='F')

x = zeros(data.GetNumberOfPoints())
y = zeros(data.GetNumberOfPoints())
z = zeros(data.GetNumberOfPoints())

for i in range(data.GetNumberOfPoints()):
        x[i],y[i],z[i] = data.GetPoint(i)

x = x.reshape(dim,order='F')
y = y.reshape(dim,order='F')
z = z.reshape(dim,order='F')