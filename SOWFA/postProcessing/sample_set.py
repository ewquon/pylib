#!/usr/bin/python
from __future__ import print_function
import sys,os
import numpy as np
from scipy.ndimage.filters import uniform_filter

def pretty_list(strlist,indent=2,sep='\t',width=80):
    """For formatting long lists of strings of arbitrary length
    """
    sep = sep.expandtabs()
    max_item_len = max([len(s) for s in strlist])
    items_per_line = int((width - (indent+max_item_len)) / (len(sep)+max_item_len)) + 1
    items_per_line = max(items_per_line,1)
    Nlines = int(len(strlist) / items_per_line)
    extraline = (len(strlist) % items_per_line) > 0
    fmtstr = '{{:{:d}s}}'.format(max_item_len)
    strlist = [ fmtstr.format(s) for s in strlist ] # pad strings so that they're all the same length
    finalline = ''
    for line in range(Nlines):
        ist = line*items_per_line
        finalline += indent*' ' + sep.join(strlist[ist:ist+items_per_line]) + '\n'
    if extraline:
        finalline += indent*' ' + sep.join(strlist[Nlines*items_per_line:]) + '\n'
    return finalline

class uniform:
    """Module for post-processing a time series of SOWFA "uniform" set samples
    (i.e., sampling along a line). Output is a simple xyz-style format.
    """

    sampleExt = 'xy'
    timeseriesfile = 'timeseries_info.npz'

    def __init__(self,path='.'):
        """Sets timeNames and sampleNames with the time directory and sample names, respectively.
        Sets t to a numpy array of times
        """
        self.path = path
        self.t = []
        self.timeNames = []
        self.sampleNames = []
        tsfile = os.path.join(path,self.timeseriesfile)

        # get list of times from directory names
        try:
            dirs = next(os.walk(path))[1]
        except StopIteration:
            print('Path',path,'not found')
            return
        times = []
        timeNames = []
        for dname in dirs:
            try:
                tval = float(dname)
                times.append( tval )
                timeNames.append( dname )
            except ValueError: pass
        self.t = np.array(times)
        self.N = len(self.t)
        if self.N==0: 
            print('No time directories found in',path)
            # attempt to load time series information if it exists
            if os.path.isfile(tsfile):
                npzfile = np.load(tsfile)
                self.t = npzfile['t']
                self.dt = npzfile['dt']
                self.timeNames = [ name.decode('utf-8') for name in npzfile['timeNames'] ]  # npzfile['timeNames']
                self.sampleNames = [ name.decode('utf-8') for name in npzfile['sampleNames'] ]  # npzfile['sampleNames']
                self.N = len(self.t)
                print('Loaded time series information from',tsfile)
            return

        # sort based on times
        order = np.argsort(self.t)
        self.t = self.t[order]
        self.timeNames = [ timeNames[i] for i in order ]
        dt = np.diff(self.t)
        if np.max(dt)-np.min(dt) > 1e-14:
            print('Warning: gaps detected in sampling times')
            self.dt = dt
        else:
            self.dt = dt[0]

        # gather sample names
        path0 = os.path.join(path,timeNames[0])
        flist = [ f for f in os.listdir(path0) if os.path.isfile(os.path.join(path0, f)) ]
        for f in flist:
            fsplit = f.split('.')
            assert( fsplit[-1] == self.sampleExt )
            name = '.'.join( fsplit[:-1] )
            self.sampleNames.append( name )

        # save time information in case we don't want to move all the unprocessed data
        np.savez(tsfile,
                 t=self.t,dt=self.dt,
                 timeNames=self.timeNames,
                 sampleNames=self.sampleNames)

    def __repr__(self):
        s = 'Read times from {:s} :\n{:s}\n'.format(self.path,str(self.t)) \
          + 'Sample names :\n' \
          + pretty_list(sorted(self.sampleNames))
        return s

    def get_sample(self,
                   name,field='U',sort=True,
                   specialnames={'p_rgh':'prgh',
                                 'ddt_halfUU':'ddtHalfUU',
                                 'div_up':'divUp',
                                 'div_utau':'divUtau'
                                },
                   oldnames=False,
                   verbose=True):
        """Returns a sampled field with the specified field name,
        assuming 'x' is identical for all samples

        For a scalar field, the output array has shape (N,NX), where NX
        is the number of spatial samples and N is the length of the time
        array; for a vector field, the output array has shape (N,NX,3).
        """
        found = False
        if field in specialnames.keys():
            searchfield = specialnames[field]
        else:
            searchfield = field

        for f in self.sampleNames:
            if not f.startswith(name): continue
            fieldstring = f[len(name):]
            for orig,new in specialnames.items():
                fieldstring = fieldstring.replace(orig,new)
            fields_in_file = fieldstring.split('_')
            if searchfield in fields_in_file:
                found = True
                pos = fields_in_file.index(searchfield)
                break
        if not found:
            print('Sample',name,'with field',field,'not found')
            return
        fname = f + '.' + self.sampleExt
        if oldnames:
            xfile = os.path.join(self.path, f+'.x.npy')
            datafile = os.path.join(self.path, f+'.'+field+'.npy')
        else:
            xfile = os.path.join(self.path, name+'.x.npy')
            datafile = os.path.join(self.path, name+'.'+field+'.npy')

        try:
            # read from pre-processed numpy data file
            x = np.load( xfile )
            data = np.load( datafile )
            self.NX = data.shape[1]
            if len(data.shape)==3:
                isVector = True
            else:
                isVector = False
            print('Data read from',datafile)

        except IOError: # read in all data (default)
            # get number of lines (i.e., number of points in sample) and sampling locations
            testfile = os.path.join(self.path, self.timeNames[0], fname)
            with open(testfile, 'r') as f:
                for i,_ in enumerate(f): pass
            self.NX = i+1

            print('Found set in {} : field {} at position {} with {} samples'.format(fname,field,pos,self.NX))

            # allocate
            x = None
            if field=='U': 
                isVector = True
                data = np.zeros( (self.N, self.NX, 3) )
            else: 
                isVector = False
                data = np.zeros( (self.N, self.NX) )
            
            # process files in all time dirs
            for it,tdir in enumerate(self.timeNames):
                sys.stderr.write('\r  reading t= %f' % self.t[it])
                data_at_time = np.loadtxt(os.path.join(self.path, tdir, fname)) 
                if x is None:
                    x = data_at_time[:,0]
                if isVector:
                    for i in range(3):
                        # column 0 is x
                        # pos="1" ==> columns 1,2,3
                        # pos="2" ==> columns 4,5,6
                        data[it,:,i] = data_at_time[:,3*(pos-1)+i+1]
                else:
                    data[it,:] = data_at_time[:,pos]
            sys.stderr.write('\n')

            print('  saving',datafile)
            try:
                np.save( xfile, x )
                np.save( datafile, data )
            except IOError as err:
                print('  warning, unable to write out npy file:',err)

        # done processing sample
        if verbose:
            print('  x min/max : [{},{}]'.format(np.min(x),np.max(x)))
            if isVector:
                print('  {} min/max : [{},{}]'.format(
                        field,np.min(data[:,:,0]),np.max(data[:,:,0])))
                print('  {} min/max : [{},{}]'.format(
                        field,np.min(data[:,:,1]),np.max(data[:,:,1])))
                print('  {} min/max : [{},{}]'.format(
                        field,np.min(data[:,:,2]),np.max(data[:,:,2])))
            else:
                print('  {} min/max : [{},{}]'.format(
                        field,np.min(data[:,:]),np.max(data[:,:])))

        # sort by x (OpenFOAM output not guaranteed to be in order...)
        if sort:
            reorder = x.argsort()
        else:
            reorder = slice(None)

        if isVector: 
            return x[reorder], data[:,reorder,:]
        else:
            return x[reorder], data[:,reorder]


class SampleCollection(object):

    def __init__(self,sampleLocations,sampledData,formatString):
        self.sampleLocations = sampleLocations  # a list of integer sampling locations for identifying file names
        self.sampledData = sampledData  # a list of sample objects
        self.formatString = formatString  # for constructing the data file name
        self.Nloc = len(sampleLocations)
        self.Nt = sampledData.N
        self.t = sampledData.t
        self.dt = sampledData.dt
        # to be filled in when sample_all is called...
        self.x = None
        self.N = None
        # to be set when calculate_means is called...
        self.tavg = None
        self.Ntavg = None

    def sample_all(self,sample_list=['U','T'],check_x=True,pointTolerance=1e-4):
        self.data = {}
        for iloc,loc in enumerate(self.sampleLocations):
            print('Sample {} at {} m'.format(iloc,loc))
            sampleName = self.formatString.format(int(loc))
            for fieldname in sample_list:
                x, dataarray = self.sampledData.get_sample(sampleName,fieldname,
                                                           sort=False,
                                                           verbose=False)
                # make sure arrays are sorted, for backwards compatibility
                #reorder = x.argsort()
                #x = x[reorder]
                #Uarray = Uarray[:,reorder,:]
                # save sampled data
                if self.x is None:
                    x.sort()
                    self.x = x
                    self.N = len(x)
                elif check_x:
                    assert(np.max(np.abs(x-self.x)) < pointTolerance)
                if fieldname not in self.data.keys():
                    #self.Udata = np.zeros((self.Nloc,self.Nt,self.N,3))
                    #self.Tdata = np.zeros((self.Nloc,self.Nt,self.N))
                    self.data[fieldname] = np.zeros([self.Nloc]+list(dataarray.shape))
                if len(dataarray.shape) == 3:
                    self.data[fieldname][iloc,:,:,:] = dataarray
                else:
                    self.data[fieldname][iloc,:,:] = dataarray

    def calculate_pressure(self,Tref,g):
        self.Tref = Tref
        self.rhok = 1.0 - ( (self.data['T'] - Tref)/Tref )  # shape == (Nloc,Nt,Nx)
        h = self.x
        print('Assuming heights are',h)
        self.p = self.p_rgh + self.rhok*g*h
        
    def calculate_means(self,tavg_window=600.0,
                        calculate_pressure=False,Tref=300.0,g=9.81):
        """
        Calculates
        ----------
        Umean, Vmean, Wmean; Tmean :
            Time-averaged velocity components, averaged over
            'tavg_window' (e.g. 10 min) [m/s] [deg K]
        ufluc, vfluc, wfluc; tfluc :
            Velocity / potential temperature fluctuations, calculated
            with umean, vmean, wmean [m/s] [deg K]
        pfluc :
            Pressure fluctuations, divided by ref density, calculated
            using the specified reference temperature Tref [m^2/s^2]
        uu, vv, ww, uv, uw, vw:
            Mean Reynolds' stress components, calculated using ufluc,
            vfluc, wfluc [m^2/s^2]
        tw:
            Heat flux, calculated using wfluc and tfluc [deg K-m/s]
        pw:
            Heat flux, calculated using wfluc and tfluc [m^3/s^3]
        _.shape = (Ntavg,N)
        """
        Navg = int(tavg_window/self.dt)
        avgrange = slice(int(Navg/2),-int(Navg/2)+1)
        self.tavg = self.t[avgrange]
        self.Ntavg = len(self.tavg) 

        #window = np.ones((Navg,))/Navg # for np.convolve

        if calculate_pressure:
            self.calculate_pressure(Tref,g)

        #
        # calculate time averages and fluctuations  
        #
        self.U_mean  = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.V_mean  = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.W_mean  = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.T_mean  = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.u_fluc  = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.v_fluc  = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.w_fluc  = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.t_fluc  = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.uu_mean = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.vv_mean = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.ww_mean = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.uv_mean = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.uw_mean = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.vw_mean = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.tw_mean = np.zeros((self.Nloc,self.Ntavg,self.N))
        self.pw_mean = np.zeros((self.Nloc,self.Ntavg,self.N))
        for iloc in range(self.Nloc):
            for ix in range(self.N):
                U = uniform_filter( self.data['U'][iloc,:,ix,0], Navg ) # size Nt
                V = uniform_filter( self.data['U'][iloc,:,ix,1], Navg )
                W = uniform_filter( self.data['U'][iloc,:,ix,2], Navg )
                T = uniform_filter( self.data['T'][iloc,:,ix], Navg )
                P = uniform_filter( self.p[iloc,:,ix], Navg )
                up = self.data['U'][iloc,:,ix,0] - U # size Nt
                vp = self.data['U'][iloc,:,ix,1] - V
                wp = self.data['U'][iloc,:,ix,2] - W
                tp = self.data['T'][iloc,:,ix] - T
                self.U_mean[iloc,:,ix] = U[avgrange] # size Ntavg
                self.V_mean[iloc,:,ix] = V[avgrange]
                self.W_mean[iloc,:,ix] = W[avgrange]
                self.T_mean[iloc,:,ix] = T[avgrange]
                self.u_fluc[iloc,:,ix] = up[avgrange]
                self.v_fluc[iloc,:,ix] = vp[avgrange]
                self.w_fluc[iloc,:,ix] = wp[avgrange]
                self.t_fluc[iloc,:,ix] = tp[avgrange]
                self.uu_mean[iloc,:,ix] = uniform_filter( up*up, Navg )[avgrange]
                self.vv_mean[iloc,:,ix] = uniform_filter( vp*vp, Navg )[avgrange]
                self.ww_mean[iloc,:,ix] = uniform_filter( wp*wp, Navg )[avgrange]
                self.uv_mean[iloc,:,ix] = uniform_filter( up*vp, Navg )[avgrange]
                self.uw_mean[iloc,:,ix] = uniform_filter( up*wp, Navg )[avgrange]
                self.vw_mean[iloc,:,ix] = uniform_filter( vp*wp, Navg )[avgrange]
                self.tw_mean[iloc,:,ix] = uniform_filter( tp*wp, Navg )[avgrange]
                self.pw_mean[iloc,:,ix] = uniform_filter( pp*wp, Navg )[avgrange]
        self.k = 0.5*(self.uu_mean + self.vv_mean + self.ww_mean)


    def calculate_TIdir(self):
        """Calculates the turbulence intensity [%] in the wind-aligned direction
        """
        Umag = np.sqrt(self.U_mean**2
                     + self.V_mean**2
                     + self.W_mean**2)
        windDir = np.abs(np.arctan2(self.V_mean,self.U_mean))
        TIdir = self.uu_mean *   np.cos(windDir)**2 \
              + self.uv_mean * 2*np.sin(windDir)*np.cos(windDir) \
              + self.vv_mean *   np.sin(windDir)**2
        self.TIdir = 100.0*TIdir/Umag

#===============================================================================
#===============================================================================
#===============================================================================
# test run
if __name__=='__main__':

    line0 = uniform(path='linesTransverse')
    print(line0)

    x,U = line0.get_sample('line_08km','U')


