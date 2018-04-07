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

if len(sys.argv) > 2:
    ny,nz = [ int(arg) for arg in sys.argv[1:3] ]
    kwargs = { 'NY': ny, 'NZ': nz }
    if len(sys.argv) > 3:
        kwargs['order'] = sys.argv[3]
    points = bc.read_boundary_points('points',**kwargs)
    y = points[0]
    z = points[1]
else:
    points = bc.read_boundary_points('points')
    y = points[0]
    z = points[1]
    ny = len(y)
    nz = len(z)

ts = TimeSeries('.')
Ntimes = len(ts.outputTimes)
print(Ntimes)
print(ts)

Uarray = np.zeros((Ntimes,ny,nz,3))
Tarray = np.zeros((Ntimes,ny,nz))
karray = np.zeros((Ntimes,ny,nz))

have_k = True
for itime,dpath in enumerate(ts.dirList):
    tname = os.path.split(dpath)[-1]
    print('t=',ts.outputTimes[itime],' ',dpath,' (',tname,')')

    Ufield = bc.read_vector_data(os.path.join(tname, 'U'),
                                 NY=ny, NZ=nz)
    for i in range(3):
        Uarray[itime,:,:,i] = Ufield[i,:,:]

    Tfield = bc.read_scalar_data(os.path.join(tname, 'T'),
                                 NY=ny, NZ=nz)
    Tarray[itime,:,:] = Tfield

    if os.path.isfile(os.path.join(tname, 'k')):
        kfield = bc.read_scalar_data(os.path.join(tname, 'k'),
                                     NY=ny, NZ=nz)
        karray[itime,:,:] = kfield
    else:
        have_k = False

npzfile = 'data.npz'
if have_k:
    np.savez_compressed(npzfile, U=Uarray, T=Tarray, k=karray)
else:
    np.savez_compressed(npzfile, U=Uarray, T=Tarray)
print('wrote',npzfile)

