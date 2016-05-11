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

#from memory_profiler import profile #-- THIS IS SLOW

class turbsim_bts:
    realtype = np.float32

    def __init__(self,prefix,Umean=None,verbose=False):
        """ Handle bts files containing binary full-field time series output from TurbSim.
        Tested with TurbSim v2.00.05c-bjj, 25-Feb-2016
        """
        self.prefix = prefix
        self.Umean = Umean
        self.hub = dict() #hub-height wind speeds
        self.field = dict() #full NY x NZ field
        self._readBTS(prefix,verbose=verbose)

    #@profile
    def _readBTS(self,prefix,verbose=False):
        """ Process AeroDyn full-field files
        """
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
            self.T  = self.realtype(self.N * self.dt)
            self.Nsize = 3*self.NY*self.NZ*self.N
            if verbose:
                print '  nt=',self.N
                print '  (problem size: {:d} points)'.format(self.Nsize)
                print '  dz,dy=',self.dz,self.dy
                print '  TimeStep=',self.dt
                print '  Period=',self.T

            # - read reference values
            self.uhub = f.read_float(dtype=self.realtype)
            self.zhub = f.read_float(dtype=self.realtype) # NOT USED
            self.zbot = f.read_float(dtype=self.realtype)
            if self.Umean is None:
                self.Umean = self.uhub
                if verbose: print '  Umean = uhub =',self.Umean,'(for calculating fluctuations)'
            else: # user-specified Umean
                if verbose:
                    print '  Umean =',self.Umean,'(for calculating fluctuations)'
                    print '  uhub=',self.uhub,' (NOT USED)'
            if verbose:
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
            # note: need to specify Fortran-order to properly read data using np.nditer
            t0 = time.clock()

            if verbose: print 'Reading normalized grid data'
            self.V = np.zeros((3,self.NY,self.NZ,self.N),order='F',dtype=self.realtype)
            print '  V size :',self.V.nbytes/1024.**2,'MB'
            for val in np.nditer(self.V, op_flags=['writeonly']):
                val[...] = f.read_int2()

            if self.Ntower > 0:
                self.Vtow = np.zeros((3,self.Ntower,self.N),order='F',dtype=self.realtype)
                print '  Vtow size :',self.Vtow.nbytes/1024.**2,'MB'
                if verbose: print 'Reading normalized tower data'
                for val in np.nditer(self.Vtow, op_flags=['writeonly']):
                    val[...] = f.read_int2()

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
            self.V[0,:,:,:] -= self.Umean

            #print '  V size :',self.V.nbytes/1042.**2,'MB'
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
            if verbose:
                #print 'Read times',self.t
                print 'Read times [',self.t[0],self.t[1],'...',self.t[-1],']'

    #--end of self._readBTS()

    #@profile
    def tileY(self,ntiles,mirror=False):
        """ Duplicate field in lateral direction
        'ntiles' is the final number of panels including the original
        Set 'mirror' to True to flip every other tile
        """
        ntiles = int(ntiles)
        print 'Creating',ntiles,'horizontal tiles'
        print '  before:',self.V.shape
        #self.V = np.tile(self.V,(1,ntiles,1,1)) #-- this repeats the ymax boundary
        if mirror:
            # [0 1 2] --> [0 1 2 1 0 1 2 .. ]
            NYnew = (self.NY-1)*ntiles + 1
            Vnew = np.zeros((3,NYnew,self.NZ,self.N))
            Vnew[:,:self.NY,:,:] = self.V[:,:self.NY,:,:]
            delta = self.NY - 1
            flipped = True
            for i in range(1,ntiles):
                if flipped:
                    Vnew[:,i*delta+1:(i+1)*delta+1,:,:] = self.V[:,delta-1::-1,:,:]
                else:
                    Vnew[:,i*delta+1:(i+1)*delta+1,:,:] = self.V[:,1:,:,:]
                flipped = not flipped
            self.V = Vnew
        else:
            # [0 1 2] --> [0 1 0 1 .. 0 1 2]
            self.V = np.tile(self.V[:,:-1,:,:],(1,ntiles,1,1))
            plane0 = np.zeros((3,1,self.NZ,self.N))
            plane0[:,0,:,:] = self.V[:,-1,:,:]
            self.V = np.concatenate((self.V,plane0),axis=1)
        print '  after :',self.V.shape


        #self.NY *= ntiles
        self.NY = NYnew
        assert( self.V.shape == (3,self.NY,self.NZ,self.N) )
        self.y = -0.5*(self.NY-1)*self.dy + np.arange(self.NY,dtype=self.realtype)*self.dy

    def writeMappedBC(self,outputdir,Uprofile=lambda z:0.0,interval=1,xinlet=0.0,Tmax=None,bcname='inlet'):
        """ For use with OpenFOAM's timeVaryingMappedFixedValue boundary condition.
        This will create a points file and time directories in 'outputdir', which should be placed in constant/boundaryData/<patchname>.
        """
        import os
        if not os.path.isdir(outputdir):
            print 'Creating output dir :',outputdir
            os.mkdir(outputdir)

        # write points file
        pointshdr = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.x                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       vectorField;
    location    "constant/boundaryData/{patchName}";
    object      points;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"""
        fname = outputdir + os.sep + 'points'
        print 'Writing',fname
        with open(fname,'w') as f:
            f.write(pointshdr.format(patchName=bcname))
            f.write('{:d}\n(\n'.format(self.NY*self.NZ))
            for j in range(self.NZ):
                for i in range(self.NY):
                    f.write('({:f} {:f} {:f})\n'.format(xinlet,self.y[i],self.z[j]))
            f.write(')\n')

        # write time dirs
        datahdr = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.x                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       vectorAverageField;
    location    "constant/boundaryData/{patchName}/{timeName}";
    object      values;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Average
(0 0 0)\n\n"""
        if Tmax is None: Tmax = self.T
        Nsteps = int(Tmax / self.dt)
        for i in range(0,Nsteps,interval):
            itime = np.mod(i,self.N)
            tname = '{:f}'.format(self.realtype(i*self.dt)).rstrip('0').rstrip('.')
            try: os.mkdir(outputdir+os.sep+tname)
            except: pass
            fname = outputdir + os.sep + tname + os.sep + 'U'
            sys.stdout.write('Writing {} (itime={})\n'.format(fname,itime))
            with open(fname,'w') as f:
                f.write(datahdr.format(patchName=bcname,timeName=tname))
                f.write('{:d}\n(\n'.format(self.NY*self.NZ))
                for j in range(self.NZ):
                    Uinlet = np.array([Uprofile(self.z[j]),0,0])
                    #print ' ',self.z[j],Uinlet
                    for i in range(self.NY):
                        #f.write('({v[0]:f} {v[1]:f} {v[2]:f})\n'.format(v=self.V[:,i,j,itime]))
                        f.write('({v[0]:f} {v[1]:f} {v[2]:f})\n'.format(v=self.V[:,i,j,itime]+Uinlet))
                f.write(')\n')


    def writeVTK(self,fname,itime=None,output_time=None,stdout='verbose'):# {{{
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
        u = np.zeros((1,self.NY,self.NZ)); u[0,:,:] = self.V[0,:,:,itime]
        v = np.zeros((1,self.NY,self.NZ)); v[0,:,:] = self.V[1,:,:,itime]
        w = np.zeros((1,self.NY,self.NZ)); w[0,:,:] = self.V[2,:,:,itime]
        VTKwriter.vtk_write_structured_points( open(fname,'wb'), #binary mode
            1,self.NY,self.NZ,
            [u,v,w],
            datatype=['vector'],
            dx=1.0,dy=self.dy,dz=self.dz,
            dataname=['fluctuations'], #dataname=['TurbSim_velocity'],
            origin=[0.,self.y[0],self.z[0]],
            indexorder='ijk')# }}}

    def writeVTKSeries(self,prefix=None,step=1,stdout='verbose'):# {{{
        """ Call writeVTK for a range of times
        """
        if not prefix: prefix = self.prefix
        for i in range(0,self.N,step):
            fname = prefix + '_' + str(i) + '.vtk'
            self.writeVTK(fname,itime=i,stdout=stdout)
	if stdout=='overwrite': sys.stdout.write('\n')# }}}

#===============================================================================
#===============================================================================
#===============================================================================
# testing
if __name__=='__main__':
    prefix = 'Kaimal_15'
    #field = turbsim_bts(prefix,verbose=True,Umean=6.8)
    field = turbsim_bts(prefix,verbose=True)

#    field.writeVTKSeries(prefix='vtk/Kaimal_15', step=10, stdout='overwrite')

#    field.writeMappedBC('west',Uprofile=lambda z:6.0,interval=10,Tmax=1000.,bcname='west')

#    field.tileY(3,mirror=True)
#    field.writeVTKSeries(prefix='vtk_tile3/Kaimal_15', step=5, stdout='overwrite')

    field.tileY(2,mirror=True)
    field.writeVTKSeries(prefix='vtk_tile2/Kaimal_15', step=10, stdout='overwrite')
    field.writeMappedBC('west2',Uprofile=lambda z:6.0,interval=10,Tmax=1200.,bcname='west')


