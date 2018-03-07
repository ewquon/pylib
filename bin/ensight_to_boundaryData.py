#!/usr/bin/env python
#
# Converts SOWFA array sampling in ensight format to boundaryData.
# Can be useful if downscaling is desired, e.g., sampling from a 20-m
# precursor to be used as inflow to a 10-m domain. Also writes out a
# compressed npz file with all of the boundary data.
#
# Note: This has only been tested for flow from the west.
#
# Written by Eliot Quon (eliot.quon@nrel.gov) -- 2018-03-07
#
import sys, os
import numpy as np
#from SOWFA.include import read_setUp
from NWTC.datatools.timeseries import TimeSeries
import NWTC.datatools.ensight as ens
import NWTC.datatools.SOWFA.timeVaryingMappedBC as bc

if len(sys.argv) <= 5:
    sys.exit('USAGE: '+sys.argv[0]+' ensight_array_outputdir BC_name flow_dir plane_x ncells_y ncells_z [new_boundaryData_dir]')
dpath = sys.argv[1]
bcname = sys.argv[2]
xinflow = sys.argv[3]
ny = int(sys.argv[4])
nz = int(sys.argv[5])
if len(sys.argv) > 6:
    boundaryData = sys.argv[6]
else:
    boundaryData = 'boundaryData'
assert(bcname == 'west')  # only handle flow from west for now
Ncells = ny * nz
Npoints = (ny+1) * (nz+1)

if os.path.isdir(boundaryData):
    print('Warning:',boundaryData,'already exists!')

#config = read_setUp()
#yscale = ny / config['ny']
#zscale = nz / config['nz']
#print('Scaling factors:',yscale,zscale)

ts = TimeSeries(dpath)
Ntimes = len(ts.outputTimes)
print(Ntimes)
print(ts)

Uarray = np.zeros((Ntimes,ny,nz,3))
Tarray = np.zeros((Ntimes,ny,nz))
karray = np.zeros((Ntimes,ny,nz))

bcdir = os.path.join(boundaryData, bcname)
if not os.path.isdir(bcdir):
    os.makedirs(bcdir)

for itime,dpath in enumerate(ts.dirList):
    print('t=',ts.outputTimes[itime])
    tname = os.path.split(dpath)[-1]
    pointsfile = os.path.join(boundaryData, bcname, 'points')
    datapath = os.path.join(boundaryData, bcname, tname)
    if not os.path.isfile(pointsfile):
        # read mesh
        # y = y.reshape((ny,nz), order='F')
        # z = z.reshape((ny,nz), order='F')
        x,y,z = ens.read_mesh(os.path.join(dpath, bcname+'_U.mesh'))
        assert(len(x) == len(y) == len(z) == Ncells)
        assert(np.min(x) == np.max(x)) # constant
        x[:] = xinflow
        bc.write_points(pointsfile, x, y, z, patchName=bcname)
        print('  wrote',pointsfile)
        np.savez_compressed(pointsfile+'.npz',
                            y=y.reshape((ny,nz), order='F'),
                            z=z.reshape((ny,nz), order='F'))
        print('  wrote',pointsfile+'.npz')
    if not os.path.isdir(datapath):
        os.makedirs(datapath)
        # process U
        u,v,w = ens.read_vector(os.path.join(dpath, bcname+'_U.000.U'), Ncells)
        bc.write_data(os.path.join(boundaryData, bcname, tname, 'U'),
                      np.stack((u,v,w)), patchName=bcname)
        Uarray[itime,:,:,0] = u.reshape((ny,nz), order='F')
        Uarray[itime,:,:,1] = v.reshape((ny,nz), order='F')
        Uarray[itime,:,:,2] = w.reshape((ny,nz), order='F')
        # process T
        T = ens.read_scalar(os.path.join(dpath, bcname+'_T_k.000.T'), Ncells)
        bc.write_data(os.path.join(boundaryData, bcname, tname, 'T'),
                      T, patchName=bcname)
        Tarray[itime,:,:] = T.reshape((ny,nz), order='F')
        # process k
        k = ens.read_scalar(os.path.join(dpath, bcname+'_T_k.000.k'), Ncells)
        bc.write_data(os.path.join(boundaryData, bcname, tname, 'k'),
                      k, patchName=bcname)
        karray[itime,:,:] = k.reshape((ny,nz), order='F')
        print('  wrote data in',datapath)
    else:
        # read existing boundaryData files at this time
        Ufield = bc.read_vector_data(os.path.join(boundaryData, bcname, tname, 'U'),
                                     NY=ny, NZ=nz)
        Tfield = bc.read_scalar_data(os.path.join(boundaryData, bcname, tname, 'T'),
                                     NY=ny, NZ=nz)
        kfield = bc.read_scalar_data(os.path.join(boundaryData, bcname, tname, 'k'),
                                     NY=ny, NZ=nz)
        for i in range(3):
            Uarray[itime,:,:,i] = Ufield[i,:,:]
        Tarray[itime,:,:] = Tfield
        karray[itime,:,:] = kfield

npzfile = os.path.join(boundaryData, bcname, bcname+'.npz')
np.savez_compressed(npzfile, U=Uarray, T=Tarray, k=karray)
print('wrote',npzfile)

