#!/usr/local/bin/python
#
# TurbSim data processing module (binary AeroDyn format)
# written by Eliot Quon (eliot.quon@nrel.gov)
#
import sys
import numpy as np
import VTKwriter
from binario import *

class turbsim_bts:
    realtype = np.float32

    def __init__(self,prefix,verbose=False):
        """ Handle bts files containing binary full-field time series output from TurbSim.
        Tested with TurbSim v2.00.05c-bjj, 25-Feb-2016
        """
        self.prefix = prefix
        self.hub = dict() #hub-height wind speeds
        self.field = dict() #full NY x NZ field
        self._readBTS(prefix,verbose=verbose)

    def _readBTS(self,prefix,verbose=False):
        fname = prefix + '.bts'
        with binaryfile(fname) as f:
            #
            # read header info
            #
            print 'Reading header information'
            ID = f.read_int2()
            assert( ID==7 or ID==8 )
            if ID==7: filetype = 'non-periodic'
            elif ID==8: filetype = 'periodic'
            else: filetype = 'UNKNOWN'
            print '  id= {:d} ({:s})'.format(ID,filetype)

            # - read resolution settings
            self.NZ = f.read_int4()
            self.NY = f.read_int4()
            self.Ntower = f.read_int4()
            if verbose:
                print '  NumGrid_Z,_Y=',self.NZ,self.NY
                print '  ntower=',self.Ntower
            self.N = f.read_int4()
            self.dz = f.read_float()
            self.dy = f.read_float()
            self.dt = f.read_float()
            if verbose:
                print '  nt=',self.N
                print '  dz,dy=',self.dz,self.dy
                print '  TimeStep=',self.dt

            # - read reference values
            self.uhub = f.read_float()
            self.zhub = f.read_float()
            self.zbot = f.read_float()
            if verbose:
                print '  uhub=',self.uhub
                print '  HubHt=',self.zhub
                print '  Zbottom=',self.zbot

            # - read scaling factors
            self.Vslope = np.zeros(3,dtype=self.realtype)
            self.Vintercept = np.zeros(3,dtype=self.realtype)
            for i in range(3):
                self.Vslope[i] = f.read_float()
                self.Vintercept[i] = f.read_float()
            if verbose:
                print '  Vslope=',self.Vslope
                print '  Vintercept=',self.Vintercept

            # - read turbsim info string
            nchar = f.read_int4()
            version = f.read(N=nchar)
            if verbose: print version

            #
            # read normalized data
            #
            #for it in range(self.N):
            #    for iz in range(self.NZ):
            #        for iy in range(self.NY):
            #            for i in range(3):
            # need to specify Fortran-order to properly read data using np.nditer
            if verbose: print 'Reading normalized grid data'
            self.V = np.zeros((3,self.NY,self.NZ,self.N),order='F',dtype=self.realtype)
            for val in np.nditer(self.V, op_flags=['writeonly']):
                val[...] = f.read_int2()

            if self.Ntower > 0:
                self.Vtow = np.zeros((3,self.Ntower,self.N),order='F',dtype=self.realtype)
                if verbose: print 'Reading normalized tower data'
                for val in np.nditer(self.Vtow, op_flags=['writeonly']):
                    val[...] = f.read_int2()
                            
            #
            # calculate dimensional velocity
            #
            if verbose: print 'Calculating velocities'
            print '  u min/max [',np.min(self.V[0,:,:,:]),np.max(self.V[0,:,:,:]),']'
            print '  v min/max [',np.min(self.V[1,:,:,:]),np.max(self.V[1,:,:,:]),']'
            print '  w min/max [',np.min(self.V[2,:,:,:]),np.max(self.V[2,:,:,:]),']'
            for i in range(3):
                self.V[i,:,:,:] -= self.Vintercept[i]
                self.V[i,:,:,:] /= self.Vslope[i]
                if self.Ntower > 0:
                    self.Vtow[i,:,:] -= self.Vintercept[i]
                    self.Vtow[i,:,:] /= self.Vslope[i]

            print '  u min/max [',np.min(self.V[0,:,:,:]),np.max(self.V[0,:,:,:]),']'
            print '  v min/max [',np.min(self.V[1,:,:,:]),np.max(self.V[1,:,:,:]),']'
            print '  w min/max [',np.min(self.V[2,:,:,:]),np.max(self.V[2,:,:,:]),']'

            #
            # calculate coordinates
            #
            if verbose: print 'Calculating coordinates'
            self.y = -0.5*(self.NY-1)*self.dy + np.arange(self.NY,dtype=self.realtype)*self.dy
            self.z = self.zbot + np.arange(self.NZ,dtype=self.realtype)*self.dz
            self.ztow = self.zbot - np.arange(self.NZ,dtype=self.realtype)*self.dz

            self.t = np.arange(self.N,dtype=self.realtype)*self.dt
            if verbose: print 'Read times',self.t

    def writeVTK(self,fname,itime=None,output_time=None,stdout='verbose'):
        """ Write out binary VTK file with a single vector field.
        Can specify time index or output time.
        Note: VTKwriter expects velocity arrays of form u[NY,NX,NZ]
        """
        if output_time:
            itime = int(output_time / self.dt)
        if not itime:
            print 'Need to specify itime or output_time'
            return
	if stdout=='overwrite':
            sys.stdout.write('\rWriting time step {:d} :  t= {:f}'.format(itime,self.t[itime]))
	else: #if stdout=='verbose':
            print 'Writing out time step',itime,': t=',self.t[itime]
        u = np.zeros((1,self.NY,self.NZ)); u[0,:,:] = self.V[0,:,:,itime]
        v = np.zeros((1,self.NY,self.NZ)); v[0,:,:] = self.V[1,:,:,itime]
        w = np.zeros((1,self.NY,self.NZ)); w[0,:,:] = self.V[2,:,:,itime]
        VTKwriter.vtk_write_structured_points( open(fname,'wb'), #binary mode
            1,self.NY,self.NZ,
            [u,v,w],
            datatype=['vector'],
            dx=1.0,dy=self.dy,dz=self.dz,
            dataname=['TurbSim_velocity'],
            origin=[0.,self.y[0],self.z[0]],
            indexorder='ijk')


    def writeVTKSeries(self,prefix=None,step=1,stdout='verbose'):
        """ Call writeVTK for a range of times
        """
        if not prefix: prefix = self.prefix
        for i in range(0,self.N,step):
            fname = prefix + '_' + str(i) + '.vtk'
            self.writeVTK(fname,itime=i,stdout=stdout)
	if stdout=='overwrite': sys.stdout.write('\n')

#===============================================================================
#===============================================================================
#===============================================================================
# testing
if __name__=='__main__':
    prefix = 'Kaimal_15'
    field = turbsim_bts(prefix,verbose=True)
    #field.writeVTKSeries(prefix='vtk/Kaimal_15') #,step=10)
    field.writeVTKSeries(prefix='vtk/Kaimal_15' ,step=10, stdout='overwrite')

