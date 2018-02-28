#!/usr/bin/env python
from __future__ import print_function
import sys, os
import numpy as np
import vtk
from vtk.util import numpy_support as VN

if len(sys.argv) <= 1:
    sys.exit('Specify VTK file(s)')

configname = 'setUp'
fieldsep = '.'

#
# find configuration file
#
configpath = None
if os.path.isfile(configname):
    configpath = configname
else:
    # search upward in dir tree
    curprefix = '../'
    while (configfile is None) and os.path.isdir(curprefix):
        if os.path.isfile(curprefix + configname):
            configpath = curprefix + configname
        else:
            curprefix += '../'
if configpath is None:
    sys.exit('Config file',configname,'not found')
else:
    print('Config file:',configpath)

#
# read grid dimensions
# (for some reason I can't seem to pull this out of the vtk file...)
#
cell_dims = np.zeros((3,),dtype=int)
with open(configpath,'r') as f:
    for line in f:
        if line.strip().startswith('nx'):
            line = line.split()
            ncells = line[1].rstrip(';')
            cell_dims[0] = int(ncells)
        elif line.strip().startswith('ny'):
            line = line.split()
            ncells = line[1].rstrip(';')
            cell_dims[1] = int(ncells)
        elif line.strip().startswith('nz'):
            line = line.split()
            ncells = line[1].rstrip(';')
            cell_dims[2] = int(ncells)
if not np.all(cell_dims > 0):
    sys.exit('Grid dimensions not properly read from',configpath,':',cell_dims)
else:
    print('Grid dimensions from',configpath,':',cell_dims)
point_dims = cell_dims[:] + 1

#
# read data
#
for fpath in sys.argv[1:]:
    if not os.path.isfile(fpath):
        print(fpath,'not found')
        continue
    dpath, prefix = os.path.split(fpath)
    prefix = os.path.splitext(prefix)[0]

    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(fpath)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    sys.stdout.write('Reading '+dpath+' '+prefix+' ...')
    sys.stdout.flush()
    reader.Update()
    print('done!')
    vdata = reader.GetOutput()

    #
    # save mesh if needed
    #
    meshpath = os.path.join(dpath, 'mesh.npz')
    if not os.path.isfile(meshpath):
        N = vdata.GetNumberOfPoints()
        assert(N == np.prod(point_dims))
        x = np.zeros(N)
        y = np.zeros(N)
        z = np.zeros(N)
        for i in range(N):
            x[i],y[i],z[i] = vdata.GetPoint(i)
        x = x.reshape(point_dims, order='F')
        y = y.reshape(point_dims, order='F')
        z = z.reshape(point_dims, order='F')
        np.savez_compressed(meshpath, x=x, y=y, z=z)
        print('- wrote mesh',meshpath)
        assert(np.all(x[:,0,0] == x[:,1,1]))
        assert(np.all(y[0,:,0] == y[1,:,1]))
        assert(np.all(z[0,0,:] == z[1,1,:]))

    #
    # process point data
    #
    pointdata = vdata.GetPointData()
    for iarray in range(pointdata.GetNumberOfArrays()):
        fieldname = pointdata.GetArrayName(iarray)
        fname = prefix + fieldsep + fieldname + '.npz'
        fieldpath = os.path.join(dpath, fname)
        varray = pointdata.GetArray(fieldname)
        Ncomp = varray.GetNumberOfComponents()
        array = VN.vtk_to_numpy(varray)
        if Ncomp == 1:
            data = array.reshape(point_dims, order='F')
            np.savez_compressed(fieldpath, data)
            print('- wrote scalar field',fieldpath)
        elif Ncomp == 3:
            shape = (point_dims[0],point_dims[1],point_dims[2],3)
            data = np.zeros(shape)
            for i in range(3):
                data[:,:,:,i] = array[:,i].reshape(point_dims, order='F')
            np.savez_compressed(fieldpath, data)
            print('- wrote vector field',fieldpath)
        else:
            shape = (point_dims[0],point_dims[1],point_dims[2],Ncomp)
            data = np.zeros(shape)
            for i in range(Ncomp):
                data[:,:,:,i] = array[:,i].reshape(point_dims, order='F')
            np.savez_compressed(fieldpath, data)
            print('- wrote field (',Ncomp,'components )',fieldpath)

