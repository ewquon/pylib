#!/usr/bin/env python
import glob
from datetime import datetime
from netCDF4 import Dataset
import numpy as np
import pickle

prefix = 'wrfout_d02_'
datestr = '2016-11-12'
met_lat, met_long = 45.638004, -120.642973  # PS-12 met mast

datafiles = glob.glob(prefix+datestr+'*')
datafiles.sort()
Nt = len(datafiles)

t = []
z = None
U = None
V = None
W = None
PB = None
PHB = None
ilat, ilong = None, None  # indices of column at the met mast location

for itime,fname in enumerate(datafiles):
    wrfdata = Dataset(fname)
    #print(fname,wrfdata.file_format)
    #print(wrfdata.variables.keys())
    #print(wrfdata['XTIME'])
    #print(wrfdata['XLAT'])
    #print(wrfdata['XLONG'])
    #print(wrfdata['PHB'])  # base-state geopotential
    #print(wrfdata['PH'])  # perturbation geopotential
    #print(wrfdata['HGT'])   # terrain height
    #print(wrfdata['U'])
    #print(wrfdata['V'])
    #print(wrfdata['W'])

    datetime_str = fname[len(prefix):]
    dt = datetime.strptime(datetime_str, '%Y-%m-%d_%H_%M_%S')
    t.append(dt)

    if ilat is None:
        # find indices of column
        xlat = wrfdata['XLAT'][:][0,:,:]  # dimension 0 is time, unlimited
        xlong = wrfdata['XLONG'][:][0,:,:]
        Nlat, Nlong = xlat.shape
        tmp = np.nonzero((xlat[:-1,:-1] <= met_lat) &
                         (xlat[1:,1:] > met_lat) &
                         (xlong[:-1,:-1] <= met_long) &
                         (xlong[1:,1:] > met_long))
        ilat = tmp[0][0]
        ilong = tmp[1][0]
        assert((met_lat >= xlat[ilat,ilong]) &
               (met_lat < xlat[ilat+1,ilong+1]) &
               (met_long >= xlong[ilat,ilong]) &
               (met_long < xlong[ilat+1,ilong+1]))
        print('Met mast at {}, {}'.format(met_lat,met_long))
        print('Selected column around [{:.4f} {:.4f}], [{:.4f} {:.4f}]'.format(xlat[ilat,ilong],
                                                                               xlat[ilat+1,ilong+1],
                                                                               xlong[ilat,ilong],
                                                                               xlong[ilat+1,ilong+1]))

    if U is None:
        Nz = wrfdata['PH'][:].shape[1] - 1
        print('Allocating arrays with size',Nt,Nz)
        U = np.zeros((Nt,Nz))
        V = np.zeros((Nt,Nz))
        W = np.zeros((Nt,Nz))
        TH = np.zeros((Nt,Nz))  # potential temperature
        PH = np.zeros((Nt,Nz+1))  # perturbation geopotential (staggered)
        PHB = np.zeros((Nt,Nz+1))  # base-state geopotential (staggered)
        z = np.zeros((Nt,Nz))  # calculated from PH, PHB, g

    # calculate height AGL of column, at cell center
    PH[itime,:] = wrfdata['PH'][:][0,:,ilat,ilong]
    PHB[itime,:] = wrfdata['PHB'][:][0,:,ilat,ilong]
    zstag = (PH[itime,:] + PHB[itime,:]) / 9.81
    assert(np.abs(zstag[0] - wrfdata['HGT'][:][0,ilat,ilong]) < 0.001)
    z[itime,:] = 0.5*(zstag[1:] + zstag[:-1]) - zstag[0]

    if itime==0:
        print('Ground geopotential, terrain height =',zstag[0],wrfdata['HGT'][:][0,ilat,ilong])
        print('Calculated column heights (AGL):',z[0,:])

    # get cell-centered velocity components
    U[itime,:] = 0.5*( wrfdata['U'][:][0,:,ilat,ilong] + wrfdata['U'][:][0,:,ilat,ilong+1] )
    V[itime,:] = 0.5*( wrfdata['V'][:][0,:,ilat,ilong] + wrfdata['V'][:][0,:,ilat+1,ilong] )
    W[itime,:] = 0.5*( wrfdata['W'][:][0,1:,ilat,ilong] + wrfdata['W'][:][0,:-1,ilat,ilong] )
    #TH[itime,:] = wrfdata['TH'][:][0,:,ilat,ilong]

    print('Processed',str(dt))
    
# save data
outfile = 'columndata_'+datestr+'.npz'
np.savez(outfile,t=t,z=z,U=U,V=V,W=W,TH=TH,PH=PH,PHB=PHB)
print('Wrote',outfile)
