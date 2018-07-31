#!/usr/bin/env python
#
# Reads in OpenFOAM boundaryData and writes out a single compressed npz
# file for all times. Run this from constant/boundaryData/BCNAME
#
# Written by Eliot Quon (eliot.quon@nrel.gov) -- 2018-03-08
#
import sys, os
import numpy as np
#from SOWFA.include import read_setUp
from NWTC.datatools.timeseries import TimeSeries
import NWTC.datatools.SOWFA.timeVaryingMappedBC as bc

npzfile = 'data.npz'

if len(sys.argv) > 2:
    ny,nz = [ int(arg) for arg in sys.argv[1:3] ]
    kwargs = { 'NY': ny, 'NZ': nz }
    if len(sys.argv) > 3:
        kwargs['order'] = sys.argv[3]
    y,z,is_structured = bc.read_boundary_points('points',**kwargs)
else:
    y,z,is_structured = bc.read_boundary_points('points')

ts = TimeSeries('.')
Ntimes = len(ts.outputTimes)
print(ts)

if is_structured:
    # regular grid
    ny = len(y)
    nz = len(z)
    N = None
    Uarray = np.zeros((Ntimes,ny,nz,3))
    Tarray = np.zeros((Ntimes,ny,nz))
    karray = np.zeros((Ntimes,ny,nz))
else:
    ny = None
    nz = None
    assert(len(y) == len(z))
    N = len(y)
    Uarray = np.zeros((Ntimes,N,3))
    Tarray = np.zeros((Ntimes,N))
    karray = np.zeros((Ntimes,N))

have_k = True
for itime,dpath in enumerate(ts.dirList):
    tname = os.path.split(dpath)[-1]
    print('t={:f} {:s} ({:s})'.format(ts.outputTimes[itime],dpath,tname))

    if is_structured:
        Ufield = bc.read_vector_data(os.path.join(tname, 'U'), NY=ny, NZ=nz)
        for i in range(3):
            Uarray[itime,:,:,i] = Ufield[i,:,:]
    else:
        Ufield = bc.read_vector_data(os.path.join(tname, 'U'))
        for i in range(3):
            Uarray[itime,:,i] = Ufield[i,:]

    if is_structured:
        Tfield = bc.read_scalar_data(os.path.join(tname, 'T'), NY=ny, NZ=nz)
        Tarray[itime,:,:] = Tfield
    else:
        Tfield = bc.read_scalar_data(os.path.join(tname, 'T'))
        Tarray[itime,:] = Tfield

    if os.path.isfile(os.path.join(tname, 'k')):
        if is_structured:
            kfield = bc.read_scalar_data(os.path.join(tname, 'k'), NY=ny, NZ=nz)
            karray[itime,:,:] = kfield
        else:
            kfield = bc.read_scalar_data(os.path.join(tname, 'k'))
            karray[itime,:] = kfield
    else:
        have_k = False

print(Uarray.shape)

if have_k:
    np.savez_compressed(npzfile, U=Uarray, T=Tarray, k=karray)
else:
    np.savez_compressed(npzfile, U=Uarray, T=Tarray)
print('wrote',npzfile)

