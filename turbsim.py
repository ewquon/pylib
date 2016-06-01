#!/usr/local/bin/python
#
# TurbSim data processing module (binary AeroDyn .bts format)
# written by Eliot Quon (eliot.quon@nrel.gov)
#
import sys
import numpy as np
import VTKwriter
from binario import *

import time

#from memory_profiler import profile #-- THIS IS SLOW
# faster to uncomment the @profile lines and then run from the command line:
#   mprof run turbsim_bts.py
#   mprof plot

class bts:
    realtype = np.float32

    def __init__(self,prefix,Umean=None,verbose=False):
        """ Handle bts files containing binary full-field time series output from TurbSim.
        Tested with TurbSim v2.00.05c-bjj, 25-Feb-2016
        """
        self.prefix = prefix
        self.Umean = Umean
        self.hub = dict() #hub-height wind speeds
        self.field = dict() #full NY x NZ field

        self.meanProfilesSet = False
        self.meanProfilesRead = False

        self.tkeProfileSet = False

        self._readBTS(prefix,verbose=verbose)

    #@profile
    def _readBTS(self,prefix,verbose=False):# {{{
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

            self.scaling = np.ones(self.NZ)

            #
            # calculate coordinates
            #
            if verbose: print 'Calculating coordinates'
            self.y = -0.5*(self.NY-1)*self.dy + np.arange(self.NY,dtype=self.realtype)*self.dy
            self.z = self.zbot + np.arange(self.NZ,dtype=self.realtype)*self.dz
            #self.ztow = self.zbot - np.arange(self.NZ,dtype=self.realtype)*self.dz #--NOT USED

            self.t = np.arange(self.N,dtype=self.realtype)*self.dt
            if verbose:
                #print 'Read times',self.t
                print 'Read times [',self.t[0],self.t[1],'...',self.t[-1],']'

    #--end of self._readBTS()# }}}

    def readMeanProfile(self,Ufile='U.dat',Vfile='V.dat',Tfile='T.dat'):# {{{
        """ Read planar averages (postprocessed separately) into arrays for interpolating
        Assumed that the heights in all files are the same
        """
        hmean,Umean,Vmean,Tmean = [], [], [], []
        with open(Ufile,'r') as f:
            for line in f:
                line = line.split()
                hmean.append( float(line[0]) )
                Umean.append( float(line[1]) )
        with open(Vfile,'r') as f:
            for line in f:
                Vmean.append( float(line.split()[1]) )
        with open(Tfile,'r') as f:
            for line in f:
                Tmean.append( float(line.split()[1]) )
        assert( len(hmean)==len(Umean) and len(Umean)==len(Vmean) and len(Vmean)==len(Tmean) )

        from scipy import interpolate
        self.Ufn = interpolate.interp1d(hmean,Umean,kind='linear',fill_value='extrapolate')
        self.Vfn = interpolate.interp1d(hmean,Vmean,kind='linear',fill_value='extrapolate')
        self.Tfn = interpolate.interp1d(hmean,Tmean,kind='linear',fill_value='extrapolate')
        self.meanProfilesRead = True

        self.setMeanProfiles( Uprofile=lambda z: [self.Ufn(z),self.Vfn(z),0.0], Tprofile=lambda z: self.Tfn(z) )
    # }}}

    def checkDivergence(self):# {{{
        print 'Calculating divergence at every point in spatio-temporal grid'
        twodx = 2.*self.uhub*self.dt
        twody = 2.*self.dy
        twodz = 2.*self.dz
        self.div = np.zeros((self.N-2,self.NY-2,self.NZ-2))
        for k in range(1,self.NZ-1):
            for j in range(1,self.NY-1):
                for i in range(1,self.N-1):
                    self.div[i-1,j-1,k-1] = \
                            (self.V[0,j  ,k  ,i+1]-self.V[0,j  ,k  ,i-1])/twodx + \
                            (self.V[1,j+1,k  ,i  ]-self.V[1,j-1,k  ,i  ])/twody + \
                            (self.V[2,j  ,k+1,i  ]-self.V[2,j  ,k-1,i  ])/twodz
        self.divmean = np.mean(field.div,axis=(1,2))
        self.divmax = np.max(np.abs(field.div),axis=(1,2))

        print 'Integrating velocity fluxes'
        dS = self.dy*self.dz
        dSy= self.dz*self.uhub*self.dt
        dSz= self.dy*self.uhub*self.dt
        self.intdiv = np.zeros(self.N-1)
        for i in range(self.N-1):
            Uin = (\
                  self.V[0, :-1, :-1,i] \
                + self.V[0, :-1,1:  ,i] \
                + self.V[0,1:  , :-1,i] \
                + self.V[0,1:  ,1:  ,i] \
                ) / 4.
            Uout = (\
                  self.V[0, :-1, :-1,i+1] \
                + self.V[0, :-1,1:  ,i+1] \
                + self.V[0,1:  , :-1,i+1] \
                + self.V[0,1:  ,1:  ,i+1] \
                ) / 4.
            Vin = (\
                  self.V[1,0, :-1,i  ] \
                + self.V[1,0, :-1,i+1] \
                + self.V[1,0,1:  ,i  ] \
                + self.V[1,0,1:  ,i+1] \
                ) / 4.
            Vout = (\
                  self.V[1,-1, :-1,i  ] \
                + self.V[1,-1, :-1,i+1] \
                + self.V[1,-1,1:  ,i  ] \
                + self.V[1,-1,1:  ,i+1] \
                ) / 4.
            Win = (\
                  self.V[2, :-1,0,i  ] \
                + self.V[2, :-1,0,i+1] \
                + self.V[2,1:  ,0,i  ] \
                + self.V[2,1:  ,0,i+1] \
                ) / 4.
            Wout = (\
                  self.V[2, :-1,-1,i  ] \
                + self.V[2, :-1,-1,i+1] \
                + self.V[2,1:  ,-1,i  ] \
                + self.V[2,1:  ,-1,i+1] \
                ) / 4.
            self.intdiv[i] = \
                dS*(np.sum(Uin) - np.sum(Uout)) + \
                dSy*(np.sum(Vin) - np.sum(Vout)) + \
                dSz*(np.sum(Win) - np.sum(Wout))# }}}

    #@profile
    def tileY(self,ntiles,mirror=False):# {{{
        """ Duplicate field in lateral direction
        'ntiles' is the final number of panels including the original
        Set 'mirror' to True to flip every other tile
        """
        ntiles = int(ntiles)
        print 'Creating',ntiles,'horizontal tiles'
        print '  before:',self.V.shape
        #self.V = np.tile(self.V,(1,ntiles,1,1)) #-- this repeats the ymax boundary
        #self.NY *= ntiles
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

        self.NY = NYnew
        assert( self.V.shape == (3,self.NY,self.NZ,self.N) )
        self.y = np.arange(self.NY,dtype=self.realtype)*self.dy
    # }}}

    def extendZ(self,zMin,zMax):# {{{
        """ Extend TurbSim domain to fit LES domain and update NZ
        Values between zMin and min(z) will be duplicated from V[:3,y,z=min(z),t]
        Values between max(z) and zMax will be zero
        """
        if zMin > self.z[0]:
            print 'zMin not changed from',self.z[0],'to',zMin
            return
        if zMax < self.z[-1]:
            print 'zMax not changed from',self.z[-1],'to',zMax
            return

        imin = int(zMin/self.dz)
        imax = int(zMax/self.dz)
        zMin = imin*self.dz
        zMax = imax*self.dz
        ioff = int((self.z[0]-zMin)/self.dz)
        print 'Extending fluctuations field in z-dir from [', self.z[0],self.z[-1],'] to [',zMin,zMax,']'
        print '  before:',self.V.shape
        
        newNZ = imax-imin+1
        Vnew = np.zeros( (3,self.NY,newNZ,self.N) )
        for iz in range(ioff):
            Vnew[:,:,iz,:] = self.V[:,:,0,:]
        Vnew[:,:,ioff:ioff+self.NZ,:] = self.V

        self.V = Vnew
        self.NZ = newNZ
        print '  after:',self.V.shape

        print 'Updating z coordinates and resetting scaling function'
        self.z = zMin + np.arange(self.NZ,dtype=self.realtype)*self.dz
        self.scaling = np.ones(self.NZ)

    # }}}

    def setMeanProfiles(self,Uprofile=lambda z:[0.0,0.0,0.0],Tprofile=lambda z:0.0):# {{{
        """ Sets the mean velocity and temperature profiles (affects output from writeMappedBC and writeVTK)
        Called by readMeanProfile after reading in post-processed planar averages.
        Can also be directly called with a user-specified analytical profile.
        """
        self.Uinlet = np.zeros((self.NZ,3))
        self.Tinlet = np.zeros(self.NZ)
        for iz,z in enumerate(self.z):
            self.Uinlet[iz,:] = Uprofile(z)
            self.Tinlet[iz]   = Tprofile(z)

        print 'Set mean profile:  z  U  T'
        for iz,U in enumerate(self.Uinlet):
            print self.z[iz],U,self.Tinlet[iz]

        self.meanProfilesSet = True
    # }}}

    def setTkeProfile(self,kprofile=lambda z:0.0):# {{{
        """ Sets the mean TKE profiles (affects output from writeMappedBC)
        Called by readHubTkeProfile.
        Can also be directly called with a user-specified analytical profile.
        """
        self.kinlet = np.zeros(self.NZ)
        for iz,z in enumerate(self.z):
            self.kinlet[iz]   = kprofile(z)

        #print 'Set TKE profile:  z  k'
        #for iz,k in enumerate(self.kinlet):
        #    print self.z[iz],k

        self.tkeProfileSet = True
    # }}}

    def setScaling(self,tanh_z90,tanh_z50,max_scaling=1.0,output=''):# {{{
        """ Set scaling of fluctuations with height
        Scaling function ranges from 0 to max_scaling. The heights at which the fluctuation magnitudes are decreased by 90% and 50% (tanh_z90 and tanh_z50, respectively) are specified to scale the hyperbolic tangent function; tanh_z90 should be set to approximately the inversion height:
          f = max_scaling * 0.5( tanh( k(z-z_50) ) + 1 )
        where
          k = arctanh(0.8) / (z_90-z_50)
        Note: If extendZ is used, that should be called to update the z coordinates prior to using this routine.
        """
        k = np.arctanh(0.8) / (tanh_z90-tanh_z50)
        self.scaling = max_scaling * 0.5*(np.tanh(-k*(self.z-tanh_z50)) + 1.0)
        fmin = np.min(self.scaling)
        fmax = np.max(self.scaling)
        assert( fmin >= 0. and fmax <= max_scaling )
        print 'Updated fluctuation scaling: range is',fmin,fmax
        
        if output:
            with open(output,'w') as f:
                f.write('# tanh scaling parameters: z_90={:f}, z_50={:f}, max_scaling={:f}\n'.format(tanh_z90,tanh_z50,max_scaling))
                f.write('# z  f(z)\n')
                for iz,z in enumerate(self.z):
                    f.write(' {:f} {:g}\n'.format(z,self.scaling[iz]))
            print 'Wrote scaling function to',output
    # }}}

    def writeMappedBC(self,outputdir='boundaryData',# {{{
            interval=1, Tstart=0., Tmax=None,
            xinlet=0.0, bcname='inlet',
            writeU=True, writeT=True, writek=True):
        """ For use with OpenFOAM's timeVaryingMappedFixedValue boundary condition.
        This will create a points file and time directories in 'outputdir', which should be placed in constant/boundaryData/<patchname>.
        The output interval is in multiples of the TurbSim time step, with output up to time Tmax.
        The inlet location and bcname should correspond to the LES inflow plane. 
        """
        import os
        if not os.path.isdir(outputdir):
            print 'Creating output dir :',outputdir
            os.makedirs(outputdir)

        if not self.meanProfilesSet:
            print 'Note: Mean profiles have not been set or read from files'
            self.setMeanProfiles()
        if writek and not self.meanProfilesSet:
            print 'Note: Mean TKE profile has not been set'
            self.setTkeProfile()

        if writeU:
            # for scaling fluctuations
            up = np.zeros((3,1,self.NY,self.NZ))

        #
        # write points file
        #
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

        #
        # write time dirs
        #
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
    class       {patchType}AverageField;
    location    "constant/boundaryData/{patchName}/{timeName}";
    object      values;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Average
{avgValue}\n\n"""
        if Tmax is None: 
            Tmax = Tstart + self.T
        istart = int(self.realtype(Tstart) / self.dt)
        iend = int(self.realtype(Tmax) / self.dt)
        print 'Outputting time length',(iend-istart)*self.dt
        for i in range(istart,iend,interval):
            itime = np.mod(i-istart,self.N)
            tname = '{:f}'.format(self.realtype(i*self.dt)).rstrip('0').rstrip('.')
            try: os.mkdir(outputdir+os.sep+tname)
            except: pass
            prefix = outputdir + os.sep + tname + os.sep

            #
            # write out U
            #
            if writeU:
                fname = prefix + 'U'
                sys.stdout.write('Writing {} (itime={})\n'.format(fname,itime))
                # scale fluctuations
                up[0,0,:,:] = self.V[0,:,:,itime]
                up[1,0,:,:] = self.V[1,:,:,itime]
                up[2,0,:,:] = self.V[2,:,:,itime]
                for iz in range(self.NZ):
                    up[:,0,:,iz] *= self.scaling[iz]
                with open(fname,'w') as f:
                    f.write(datahdr.format(patchType='vector',patchName=bcname,timeName=tname,avgValue='(0 0 0)'))
                    f.write('{:d}\n(\n'.format(self.NY*self.NZ))
                    for k in range(self.NZ):
                        for j in range(self.NY):
                            #f.write('({v[0]:f} {v[1]:f} {v[2]:f})\n'.format(v=self.V[:,j,k,itime]+self.Uinlet[k,:]))
                            f.write('({v[0]:f} {v[1]:f} {v[2]:f})\n'.format(v=self.Uinlet[k,:]+up[:,0,j,k]))
                    f.write(')\n')

            #
            # write out T
            #
            if writeT:
                fname = prefix + 'T'
                sys.stdout.write('Writing {} (itime={})\n'.format(fname,itime))
                with open(fname,'w') as f:
                    f.write(datahdr.format(patchType='scalar',patchName=bcname,timeName=tname,avgValue='0'))
                    f.write('{:d}\n(\n'.format(self.NY*self.NZ))
                    for j in range(self.NZ):
                        for i in range(self.NY):
                            f.write('{s:f}\n'.format(s=self.Tinlet[j]))
                    f.write(')\n')

            #
            # write out k
            #
            if writek:
                fname = prefix + 'k'
                sys.stdout.write('Writing {} (itime={})\n'.format(fname,itime))
                with open(fname,'w') as f:
                    f.write(datahdr.format(patchType='scalar',patchName=bcname,timeName=tname,avgValue='0'))
                    f.write('{:d}\n(\n'.format(self.NY*self.NZ))
                    for j in range(self.NZ):
                        for i in range(self.NY):
                            f.write('{s:f}\n'.format(s=self.kinlet[j]))
                    f.write(')\n')
    # }}}

    def writeVTK(self,fname,itime=None,output_time=None,stdout='verbose'):# {{{
        """ Write out binary VTK file with a single vector field.
        Can specify time index or output time.
        """
        if not self.meanProfilesSet: self.setMeanProfiles()

        if output_time:
            itime = int(output_time / self.dt)
        if itime is None:
            print 'Need to specify itime or output_time'
            return
	if stdout=='overwrite':
            sys.stdout.write('\rWriting time step {:d} :  t= {:f}'.format(itime,self.t[itime]))
	else: #if stdout=='verbose':
            print 'Writing out VTK for time step',itime,': t=',self.t[itime]

        # scale fluctuations
        up = np.zeros((1,self.NY,self.NZ)); up[0,:,:] = self.V[0,:,:,itime]
        vp = np.zeros((1,self.NY,self.NZ)); vp[0,:,:] = self.V[1,:,:,itime]
        wp = np.zeros((1,self.NY,self.NZ)); wp[0,:,:] = self.V[2,:,:,itime]
        for iz in range(self.NZ):
            up[0,:,iz] *= self.scaling[iz]
            vp[0,:,iz] *= self.scaling[iz]
            wp[0,:,iz] *= self.scaling[iz]

        # calculate instantaneous velocity
        U = up.copy()
        V = vp.copy()
        W = wp.copy()
        for iz in range(self.NZ):
            U[0,:,iz] += self.Uinlet[iz,0]
            V[0,:,iz] += self.Uinlet[iz,1]
            W[0,:,iz] += self.Uinlet[iz,2]

        # write out VTK
        VTKwriter.vtk_write_structured_points( open(fname,'wb'), #binary mode
            1,self.NY,self.NZ,
            [ U,V,W, up,vp,wp ],
            datatype=['vector','vector'],
            dx=1.0,dy=self.dy,dz=self.dz,
            dataname=['U','u\''], #['fluctuations'], #dataname=['TurbSim_velocity'],
            origin=[0.,self.y[0],self.z[0]],
            indexorder='ijk')
        # }}}

    def writeVTKSeries(self,outputdir='.',prefix=None,step=1,stdout='verbose'):# {{{
        """ Driver for writeVTK to output a range of times
        """
        if not prefix: prefix = self.prefix
        import os
        if not os.path.isdir(outputdir):
            print 'Creating output dir :',outputdir
            os.makedirs(outputdir)

        for i in range(0,self.N,step):
            fname = outputdir + os.sep + prefix + '_' + str(i) + '.vtk'
            self.writeVTK(fname,itime=i,stdout=stdout)
	if stdout=='overwrite': sys.stdout.write('\n')
    # }}}

#===============================================================================
#===============================================================================
#===============================================================================
# testing
if __name__=='__main__':
    prefix = 'Kaimal_15'
    #field = turbsim_bts(prefix,verbose=True,Umean=6.8)
    field = turbsim_bts(prefix,verbose=True)

    field.checkDivergence()

#    field.writeVTKSeries(prefix='vtk/Kaimal_15', step=10, stdout='overwrite')
#    field.writeMappedBC('west',Uprofile=lambda z:6.0,interval=10,Tmax=1000.,bcname='west')

#    field.tileY(3,mirror=True)
#    field.writeVTKSeries(prefix='vtk_tile3/Kaimal_15', step=5, stdout='overwrite')

#    field.tileY(2,mirror=True)
#    field.writeVTKSeries(prefix='vtk_tile2/Kaimal_15', step=10, stdout='overwrite')
#    field.writeMappedBC('west2',Uprofile=lambda z:6.0,interval=10,Tmax=1200.,bcname='west')


