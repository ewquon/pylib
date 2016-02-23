#!/usr/bin/python
import sys
import os
import glob
import shlex #for smart splitting to preserve enquoted text
import numpy as np
import matplotlib.pyplot as plt

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

    def __init__(self,searchDir='.',searchStr='*.csv',dt=1.0):
        self.index      = 0
        self.searchDir  = searchDir
        self.searchStr  = searchStr
        self.dt         = dt
        self.N          = -1
        self.times      = []
        self.filelist   = []
        self.latest     = ''
        self.data       = dict()

        self._find_files()

    def _find_files(self):
        '''Search through searchDir for csv files
        '''
        flist = glob.glob(self.searchDir+os.sep+self.searchStr)
        self.N = len(flist)
        self.times = np.zeros((self.N))
        for i,csvfile in enumerate(flist):
            t = float( csvfile[:-4].split('_')[-1] ) * self.dt
            self.times[i] = t
        reorder = [ i[0] for i in sorted(enumerate(self.times), key=lambda x:x[1]) ]
        self.times = self.times[reorder]
        self.filelist = [ flist[i] for i in reorder ]
        self.latest = self.filelist[-1]

    def __repr__(self):
        '''Custom output
        '''
        output = 'Series containing '+str(self.N)+' data files in \''+self.searchDir+'\':\n' + \
                '  time range {:f} to {:f}  (dt={:f})'.format(self.times[0],self.times[-1],self.dt)
        if verbose:
            for i,t in enumerate(self.times):
                output += '\n  t={:f} : {:s}'.format(t,self.filelist[i])
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
        print 'Processing all files in',self.searchDir
        for i,fname in enumerate(self.filelist):
            self.data[i] = csvfile(fname,**kwargs)
            if verbose: print 'Processed',self.data[i]

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

