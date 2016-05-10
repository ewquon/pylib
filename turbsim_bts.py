#!/usr/local/bin/python
#
# TurbSim data processing module (binary AeroDyn format)
# written by Eliot Quon (eliot.quon@nrel.gov)
#
import sys
import numpy as np
import VTKwriter
from binario import *

import time

# THIS SLOWS DOWN THE FILE I/O BY AT LEAST A FACTOR OF 2x
#try:
#    import progressbar
#    showprogress = True
#except:
#    showprogress = False
showprogress = False


class turbsim_bts:
    realtype = np.float32

    def __init__(self,prefix,verbose=False,Umean=0.0):
        """ Handle bts files containing binary full-field time series output from TurbSim.
        Tested with TurbSim v2.00.05c-bjj, 25-Feb-2016
        """
        self.prefix = prefix
        self.Umean = Umean
        self.hub = dict() #hub-height wind speeds
        self.field = dict() #full NY x NZ field
        self._readBTS(prefix,verbose=verbose)

    def _readBTS(self,prefix,verbose=False):
        """ Process AeroDyn full-field files
        """
        if verbose: print 'Umean =',self.Umean
        fname = prefix + '.bts'
        with binaryfile(fname) as f:
            #
            # read header info
            #
            if verbose: print 'Reading header information from',fname

            ID = f.read_int2()
            assert( ID==7 or ID==8 )
            if ID==7: filetype = 'non-periodic'
            elif ID==8: filetype = 'periodic'
            else: filetype = 'UNKNOWN'
            if verbose: print '  id= {:d} ({:s})'.format(ID,filetype)

            # - read resolution settings
            self.NZ = f.read_int4()
            self.NY = f.read_int4()
            self.Ntower = f.read_int4()
            if verbose:
                print '  NumGrid_Z,_Y=',self.NZ,self.NY
                print '  ntower=',self.Ntower
            self.N = f.read_int4()
            self.dz = f.read_float(dtype=self.realtype)
            self.dy = f.read_float(dtype=self.realtype)
            self.dt = f.read_float(dtype=self.realtype)
            self.Nsize = 3*self.NY*self.NZ*self.N
            if verbose:
                print '  nt=',self.N
                print '  (problem size: {:d})'.format(self.Nsize)
                print '  dz,dy=',self.dz,self.dy
                print '  TimeStep=',self.dt

            # - read reference values
            self.uhub = f.read_float(dtype=self.realtype) # NOT USED
            self.zhub = f.read_float(dtype=self.realtype) # NOT USED
            self.zbot = f.read_float(dtype=self.realtype)
            if verbose:
                print '  uhub=',self.uhub,' (NOT USED)'
                print '  HubHt=',self.zhub,' (NOT USED)'
                print '  Zbottom=',self.zbot

            # - read scaling factors
            self.Vslope = np.zeros(3,dtype=self.realtype)
            self.Vintercept = np.zeros(3,dtype=self.realtype)
            for i in range(3):
                self.Vslope[i] = f.read_float(dtype=self.realtype)
                self.Vintercept[i] = f.read_float(dtype=self.realtype)
            if verbose:
                # output is float64 precision by default...
                #print '  Vslope=',self.Vslope
                #print '  Vintercept=',self.Vintercept
                print '  Vslope=',[self.Vslope[i] for i in range(3)]
                print '  Vintercept=',[self.Vintercept[i] for i in range(3)]

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
            t0 = time.clock()

            if verbose:
                if showprogress:
                    pbardesc = ['Reading normalized grid data ',progressbar.Percentage(),progressbar.Bar()]
                    pbar = progressbar.ProgressBar( widgets=pbardesc, maxval=self.Nsize).start()
                else: print 'Reading normalized grid data'
            self.V = np.zeros((3,self.NY,self.NZ,self.N),order='F',dtype=self.realtype)
            for ival,val in zip( np.arange(self.Nsize), np.nditer(self.V, op_flags=['writeonly']) ):
                val[...] = f.read_int2()
                if verbose and showprogress: pbar.update(ival)
            if verbose and showprogress: pbar.finish()

            if self.Ntower > 0:
                self.Vtow = np.zeros((3,self.Ntower,self.N),order='F',dtype=self.realtype)
                ival = 0
                if verbose:
                    if showprogress:
                        pbardesc = ['Reading normalized tower data ',progressbar.Percentage(),progressbar.Bar()]

                        pbar = progressbar.ProgressBar( widgets=pbardesc, maxval=3*self.Ntower*self.N).start()
                    else: print 'Reading normalized tower data'
                for val in np.nditer(self.Vtow, op_flags=['writeonly']):
                    val[...] = f.read_int2()
                    if verbose and showprogress: pbar.update(ival)
                    ival += 1
                if verbose and showprogress: pbar.finish()

            if verbose: print '  Read velocitiy fields in',time.clock()-t0,'s'
                            
            #
            # calculate dimensional velocity
            #
            if verbose: print 'Calculating velocities'
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
    #--end of self._readBTS

    def tileY(self,ntiles):
        """ Duplicate field in lateral direction
        ntiles is the final number of panels including the original
        """
        ntiles = int(ntiles)
        print 'Creating',ntiles,'horizontal tiles'
        print '  before:',self.V.shape
        self.V = np.tile(self.V,(1,ntiles,1,1))
        print '  after :',self.V.shape
        self.NY *= ntiles
        assert( self.V.shape == (3,self.NY,self.NZ,self.N) )

    def writeVTK(self,fname,itime=None,output_time=None,stdout='verbose'):
        """ Write out binary VTK file with a single vector field.
        Can specify time index or output time.
        """
        if output_time:
            itime = int(output_time / self.dt)
        if itime is None:
            print 'Need to specify itime or output_time'
            return
	if stdout=='overwrite':
            sys.stdout.write('\rWriting time step {:d} :  t= {:f}'.format(itime,self.t[itime]))
	else: #if stdout=='verbose':
            print 'Writing out time step',itime,': t=',self.t[itime]
        u = np.zeros((1,self.NY,self.NZ)); u[0,:,:] = self.V[0,:,:,itime] - self.Umean
        v = np.zeros((1,self.NY,self.NZ)); v[0,:,:] = self.V[1,:,:,itime]
        w = np.zeros((1,self.NY,self.NZ)); w[0,:,:] = self.V[2,:,:,itime]
        VTKwriter.vtk_write_structured_points( open(fname,'wb'), #binary mode
            1,self.NY,self.NZ,
            [u,v,w],
            datatype=['vector'],
            dx=1.0,dy=self.dy,dz=self.dz,
            dataname=['fluctuations'], #dataname=['TurbSim_velocity'],
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
    field = turbsim_bts(prefix,verbose=True,Umean=6.8)

    field.tileY(3)

    #field.writeVTKSeries(prefix='vtk/Kaimal_15') #,step=10)
    #field.writeVTKSeries(prefix='vtk/Kaimal_15', step=5, stdout='overwrite')
    field.writeVTKSeries(prefix='vtk_tile3/Kaimal_15', step=5, stdout='overwrite')

