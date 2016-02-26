#!/usr/bin/python
import sys
import os
import glob
import shlex #for smart splitting to preserve enquoted text
import numpy as np
import matplotlib.pyplot as plt
import pickle #for storing previous processed data

verbose = False

class csvfile:
    '''Data structure to hold data from csv file
    '''

    def __init__(self,fname,varnames=[],sort=False,removeDuplicates=False):
        self.fname      = fname
        self.varnames   = varnames
        self.Noutputs   = -1
        self.N          = -1

        self.x = []
        self.y = []
        self._read()
        if sort: self._sort()
        if removeDuplicates: self._remove_dup()

    def __repr__(self):
        '''Custom output
        '''
        if verbose:
            output = '{:s} : {:d} data points with variables:'.format(self.fname,self.N)
            for varname in self.varnames:
                output += '\n  '+varname
        else:
            output = '{:s} : {:d} data points'.format(self.fname,self.N)
        output += '\n'
        return output

    def _read(self):
        # first count the number of data points so we can pre-allocate
        #   and process the header if necessary
        with open(self.fname,'r') as f:
            #hdr = f.readline().split(',')
            hdrline = shlex.shlex( f.readline() )
            hdrline.whitespace += ','
            hdr = list(hdrline)
            for i,line in enumerate(f): pass
        N = i+1 # number of data points in the file
        Noutputs = len(hdr) - 1 # don't count the independent variable
        self.N = N
        self.Noutputs = Noutputs
        if len(self.varnames)==0:
            self.varnames = [ s.replace('"','') for s in hdr ]
        else: assert( len(self.varnames)==Noutputs )

        # then read the data
        self.x = np.zeros((N))
        self.y = np.zeros((N,Noutputs))
        with open(self.fname,'r') as f:
            f.readline() # skip header
            for i,line in enumerate(f):
                data = [ float(val) for val in line.split(',') ]
                self.x[i] = data[0]
                self.y[i,:] = data[1:]

    def _sort(self):
        if verbose: print '- sorting'
        order = self.x.argsort()
        self.x = self.x[order]
        for j in range(self.Noutputs):
            self.y[:,j] = self.y[order,j]

    def _remove_dup(self):
        if verbose: print '- checking for duplicates'
        dx = np.diff(self.x)
        check = np.abs(dx) < 1e-6
        if any(check):
            print 'Warning:',np.count_nonzero(check),'duplicate(s) found in',self.fname
            N = self.N
            x = self.x
            y = self.y
            newx = np.zeros((N))
            newy = np.zeros((N,self.Noutputs))
            norm = np.ones((N))
            newN = N
            inew = 0
            newx[0] = x[0]
            newy[0,:] = y[0,:]
            for i in range(1,N):
                if x[i] - newx[inew] < 1e-6: #duplicate x
                    if verbose: print 'DUPLICATE x ({:f},{:f}) == ({:f},{:f})'.format(x[inew],y[inew,0],x[i],y[i,0])
                    newN -= 1
                    newy[inew,:] += y[i,:]
                    norm[inew] += 1
                else:
                    inew += 1
                    newx[inew] = x[i]
                    newy[inew,:] = y[i,:]
            if verbose:
                print '  number of duplicated points :',np.count_nonzero(norm>1)
                print '  maximum number overlapping  :',np.max(norm)
            self.x = newx[:newN]
            for j in range(self.Noutputs):
                newy[:newN,j] /= norm[:newN]
            self.y = newy[:newN,:]
            #for xi,yi in zip(self.x,self.y[:,0]): print xi,yi #DEBUG

#-------------------------------------------------------------------------------
class series:
    '''Data structure to find and store information about a data series
    '''

    def __init__(self,searchDir='.',searchStr='*.csv',dt=1.0,histGapAction=''):
        self.index      = 0
        self.searchDir  = searchDir
        self.searchStr  = searchStr
        self.dt         = dt
        self.N          = -1
        self.times      = []
        self.filelist   = []
        self.latest     = ''
        self.data       = dict()

        self.histGapAction = histGapAction

        self._find_files()

    def _find_files(self):
        '''Search through searchDir for csv files
        '''
        flist = glob.glob(self.searchDir+os.sep+self.searchStr)
        self.N = len(flist)
        self.times = np.zeros((self.N))
        self.indices = np.zeros((self.N),dtype=int)
        for i,csvfile in enumerate(flist):
            idx = csvfile[:-4].split('_')[-1]
            t = float( idx ) * self.dt
            self.indices[i] = idx
            self.times[i] = t
        reorder = [ i[0] for i in sorted(enumerate(self.times), key=lambda x:x[1]) ]

        # checking for gaps in data
        # assuming the sampling interval is a constant
        self.indices = self.indices[reorder]
        idx_delta = self.indices[1] - self.indices[0]
        gaps = np.nonzero( np.diff(self.indices) > idx_delta )[0]
        imax = None
        if len(gaps) > 0:
            print 'Note:',len(gaps),'gaps in output history found'
            # handle cases with missing data
            if self.histGapAction=='truncate':
                imax = gaps[0]+1
                print '      truncating series at t >= {:f} (N reduced from {:d} to {:d})'.format(
                        self.times[reorder[imax]], self.N, imax )
                self.N = imax

        self.times = self.times[reorder[:imax]]
        self.filelist = [ flist[i] for i in reorder[:imax] ]
        self.latest = self.filelist[-1]

    def __repr__(self):
        '''Custom output
        '''
        output = 'Series containing '+str(self.N)+' data files in \''+self.searchDir+'\':\n' + \
                '  time range {:f} to {:f}  (dt={:f})'.format(self.times[0],self.times[-1],self.dt)
        if verbose:
            for i,t in enumerate(self.times):
                output += '\n  t={:f} : {:s}'.format(t,self.filelist[i])
        output += '\n'
        return output

    def __iter__(self):
        '''Provide iteration functionality
        '''
        return self
    def next(self):
        try:
            t = self.times[self.index]
            f = self.filelist[self.index]
            self.index += 1
        except IndexError:
            self.index = 0
            raise StopIteration
        return t,f

    def process_all(self,**kwargs):
        '''Process all files in search directory; keyword arguments are passed to the csv reader for each file
        '''
        try:
            existing_filelist = pickle.load( open('waterSurface_filelist.pkl','r') )
            print 'Found previous file list'
            filesAreTheSame = True
            if len(existing_filelist) == len(self.filelist):
                for f1,f2 in zip(existing_filelist, self.filelist):
                    if not f1==f2:
                        filesAreTheSame = False
                        break
                if filesAreTheSame: 
                    print 'Files have not changed, loading pickled data'
                    try: 
                        self.data = pickle.load( open('waterSurface_data.pkl','r') )
                        return
                    except (IOError,EOFError):
                        print 'Problem loading waterSurface_data.pkl'
                else:
                    print 'File list has changed'
        except (IOError,EOFError): pass

        pickle.dump( self.filelist, open('waterSurface_filelist.pkl','w') )

        print 'Processing all files in',self.searchDir
        for i,fname in enumerate(self.filelist):
            self.data[i] = csvfile(fname,**kwargs)
            if verbose: print 'Processed',self.data[i]

        pickle.dump( self.data, open('waterSurface_data.pkl','w') )


    def sample(self,xval=0.0,yvar=0,method='linear',plot=False):
        '''Sample a variable at a location over time, assuming all files have the same number of variables and the data contained within is sorted.
        Returns time and data arrays
        '''
        print 'Extracting',self.data[0].varnames[yvar],'at x=',xval,'using',method,'method'
        sample = np.zeros((self.N))
        #Nexact = 0
        for i in range(self.N):

            if method=='nearest':
                delta = np.abs( xval - self.data[i].x )
                imin = np.argmin(delta)
                sample[i] = self.data[i].y[imin,yvar]
                #if delta[imin]==0: Nexact += 1

            elif method=='linear': #assume data is sorted
                delta = xval - self.data[i].x
                if np.any(delta==0):
                    ival = np.nonzero(delta==0)[0][0]
                    sample[i] = self.data[i].y[ival,yvar]
                    #Nexact += 1
                else:
                    sgnchange = np.sign(delta[:-1])*np.sign(delta[1:])
                    idx = np.nonzero(delta<0)[0][0]
                    y0 = self.data[i].y[idx,yvar]
                    y1 = self.data[i].y[idx+1,yvar]
                    x0 = self.data[i].x[idx]
                    x1 = self.data[i].x[idx+1]
                    sample[i] = y0 + (y1-y0)/(x1-x0)*delta[idx]

            else:
                print 'Unrecognized method'
                break

        #print ' ',Nexact,'/',self.N,'exact matches'

        if plot:
            #plt.figure()
            plt.plot(self.times,sample,linewidth=2,label='x='+str(xval))
            plt.xlabel('time (s)')
            plt.ylabel(self.data[0].varnames[yvar])
            if method=='nearest': plt.suptitle('sampled near x='+str(xval)+' m')
            else: plt.suptitle('sampled at x='+str(xval)+' m')
            #plt.show()

        return self.times, sample

#-------------------------------------------------------------------------------
class spectrum:
    '''Data structure to store spectral components
    '''

    def __init__(self,coefficients):
        self.fname = coefficients
        self.N = -1     # number of wave components
        self.S = None   # spectral amplitude, m^2/(rad/s)
        self.w = None   # frequency, rad/s
        self.dw = -1    # frequency separation of each component, rad/s
        self.k = None   # wavenumber, 1/m
        self.p = None   # phases, rad
        self.A = None   # component amplitude, m
        self.T = None   # component period, s
        self.L = None   # component wavelength, m
        self.header = dict()

        self._read()

    def _read(self,verbose=False):
        print 'Reading spectral coefficients from',self.fname
        #--first pass
        try:
            with open(self.fname,'r') as f:
                if self.N < 0:
                    Nhdr = 0
                    for i,line in enumerate(f):
                        if line.startswith('#'): Nhdr += 1
                    self.N = i+1-Nhdr
        except IOError:
            print 'problem opening file!'
            return
        self.w = np.zeros((self.N))
        self.S = np.zeros((self.N))
        self.p = np.zeros((self.N))
        self.k = np.zeros((self.N))
        #--second pass
        iwave = 0
        with open(self.fname,'r') as f:
            for line in f:
                if line.strip()=='': break # stop on blank line
                if line.startswith('#'): 
                    if verbose: sys.stdout.write(line)
                    # parse additional parameters...
                    line = line.split()
                    try:
                        if not line[1][-1] == ':': continue
                        param = line[1][:-1]
                        self.header[param] = float(line[2])
                    except IndexError: pass
                    continue
                # z(x,t) = dw * np.sum( A*np.cos( k*x - w*(t-toffset) - phi ) )
                w, S, phi, k = [ float(val) for val in line.split() ]
                self.w[iwave] = w   
                self.S[iwave] = S
                self.p[iwave] = phi   
                self.k[iwave] = k
                iwave += 1
        self.T = 2.0*np.pi/w
        self.L = 2.0*np.pi/k
        self.A = self.S**0.5

    def __repr__(self):
        '''Custom output
        '''
        output = 'Spectrum containing {:d} wave components\n'.format(self.N)
        output+= 'Read from file: {:s}\n'.format(self.fname)
        output+= 'Header data:\n'
        for param in self.header.keys():
            output+= '  {:s}: {:f}\n'.format(param,self.header[param])
        return output


