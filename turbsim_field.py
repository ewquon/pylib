#!/usr/local/bin/python
#
# TurbSim data processing module
# written by Eliot Quon (eliot.quon@nrel.gov)
#
import numpy as np
import VTKwriter


class turbsim_field:
    def __init__(self,prefix,verbose=False):
        """ Handle set of files <prefix>.u, <prefix>.v, and <prefix>.w containing a full-field formatted time series output from TurbSim.
        Tested with TurbSim v2.00.05c-bjj, 25-Feb-2016
        """
        self.prefix = prefix
        self.hub = dict() #hub-height wind speeds
        self.field = dict() #full NY x NZ field
        self._readTurbSimScalar(prefix,'u',verbose=verbose)
        self._readTurbSimScalar(prefix,'v',verbose=verbose)
        self._readTurbSimScalar(prefix,'w',verbose=verbose)

    def _readTurbSimScalar(self,prefix,varname,verbose=False):
        fname = prefix + '.' + varname
        if verbose: print '[',fname,']'
        with open(fname,'r') as f:
            #
            # read header info
            #
            f.readline() #blank
            line = f.readline() #note
            if verbose: print line.rstrip()
            f.readline() #blank
            line = f.readline() #header
            if verbose: print line.rstrip()
            line = f.readline() #parameters
            if verbose: print line.rstrip()
            params = line.split()
            NY = int(params[0])
            NZ = int(params[1])
            dy = float(params[2])
            dz = float(params[3])
            dt = float(params[4])
            zhub = float(params[5])
            Uinf = float(params[6])
            try: assert( self.NY == NY )
            except AttributeError: self.NY = NY
            try: assert( self.NZ == NZ )
            except AttributeError: self.NZ = NZ
            try: assert( self.dy == dy )
            except AttributeError: self.dy = dy
            try: assert( self.dz == dz )
            except AttributeError: self.dz = dz
            try: assert( self.dt == dt )
            except AttributeError: self.dt = dt
            try: assert( self.zhub == zhub )
            except AttributeError: self.zhub = zhub
            try: assert( self.Uinf == Uinf )
            except AttributeError: self.Uinf = Uinf
            #
            # read coordinates
            #
            f.readline() #blank
            f.readline() #z-coordinates [m]
            z = [ float(val) for val in f.readline().split() ]
            f.readline() #blank
            f.readline() #y-coordinates [m]
            y = [ float(val) for val in f.readline().split() ]
            try: assert( self.y == y )
            except AttributeError: self.y = y
            try: assert( self.z == z )
            except AttributeError: self.z = z
            # 
            # count time steps
            #
            for i,line in enumerate(f): pass
            nlines = i+1
            assert( np.mod(nlines,self.NZ+2)==0 )
            N = nlines / (self.NZ+2)
            try: assert( self.N == N )
            except AttributeError: 
                self.N = N
                self.t = np.zeros(N)
            self.hub[varname] = np.zeros(N)
            self.field[varname] = np.zeros((N,NZ,NY))
        with open(fname,'r') as f:
            for i in range(11): f.readline()
            #
            # read time steps
            #
            for itime in range(N):
                f.readline() #blank
                line = f.readline().split() #time, hub-height speed
                self.t[itime] = float(line[0])
                self.hub[varname][itime] = float(line[1])
                for iz in range(self.NZ):
                    self.field[varname][itime,iz,:] = [ float(val) for val in f.readline().split() ]
            if verbose: print 'Read times',self.t


    def writeVTK(self,fname,itime=None,output_time=None):
        """ Write out binary VTK file with a single vector field.
        Can specify time index or output time.
        Note: VTKwriter expects velocity arrays of form u[NY,NX,NZ]
        """
        if output_time:
            itime = int(output_time / self.dt)
        if not itime:
            print 'Need to specify itime or output_time'
            return
        print 'Writing out time step',itime,': t=',self.t[itime]
        u = np.zeros((self.NY,1,self.NZ)); u[:,0,:] = np.flipud(self.field['u'][itime,:,:]).T
        v = np.zeros((self.NY,1,self.NZ)); v[:,0,:] = np.flipud(self.field['v'][itime,:,:]).T
        w = np.zeros((self.NY,1,self.NZ)); w[:,0,:] = np.flipud(self.field['w'][itime,:,:]).T
        VTKwriter.vtk_write_structured_points( open(fname,'wb'), #binary mode
            1,self.NY,self.NZ,
            [u,v,w],
            datatype=['vector'],
            dx=1.0,dy=self.dy,dz=self.dz,
            dataname=['TurbSim_velocity'],
            origin=[0.,self.y[0],self.z[0]] )


    def writeVTKSeries(self,prefix=None,step=1):
        """ Call writeVTK for a range of times
        """
        if not prefix: prefix = self.prefix
        for i in range(0,self.N,step):
            fname = prefix + '_' + str(i) + '.vtk'
            self.writeVTK(fname,itime=i)

#===============================================================================
#===============================================================================
#===============================================================================
# testing
if __name__=='__main__':
    prefix = 'Kaimal_15'
    field = turbsim_field(prefix,verbose=True)
    field.writeVTKSeries()
