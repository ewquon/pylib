#!/usr/bin/env python
#
# Helper module for postprocessing FAST outputs
# Tested with OpenFAST v8.17.00a-bjj
#
# Written by Eliot Quon (eliot.quon@nrel.gov) 2017-07-31
#
import numpy as np
import matplotlib.pyplot as plt

def read(outputFile,**kwargs):
    return FASToutput(outputFile,**kwargs)

class FASToutput(object):
    def __init__(self,fname=None,verbose=True):
        # inputs
        self.fname = fname
        self.verbose = verbose
        # initialize members
        self.outputs = None
        self.units = None
        self.Noutputs = 0
        self.N = 0
        # read output file
        if fname is not None:
            self._readFASToutput(fname)
            if verbose:
                self.printStats()

    def _readFASToutput(self,fname,Nheaderlines=6):
        if self.verbose: print 'Reading header info from',fname
        with open(fname,'r') as f:
            if self.verbose:
                for _ in range(Nheaderlines): print f.readline().strip()
            else:
                for _ in range(Nheaderlines): f.readline()
            # read names of output quantities
            self.outputs = f.readline().split()
            # read units for each output quantity
            self.units = [ s.strip('()') for s in f.readline().split() ]
            self.Noutputs = len(self.outputs)
            assert(self.Noutputs == len(self.units))
            for iline,_ in enumerate(f): pass
        self.N = iline + 1
        # read data
        if self.verbose: print 'Reading data...'
        data = np.loadtxt(fname,skiprows=Nheaderlines+2)
        for i,output in enumerate(self.outputs):
            setattr(self,output,data[:,i])
        assert(len(self.Time) == self.N)
        # aliases
        self._setAlias('t','Time')
        self._setAlias('P','RotPwr')
        self._setAlias('T','RotThrust','LSShftFxa','LSShftFxs','LSSGagFxa','LSSGagFxs')

    def _setAlias(self,name,*aliases):
        for alias in aliases:
            try:
                getattr(self, alias)
            except AttributeError:
                continue
            setattr(self, name, getattr(self,alias))
            if self.verbose: print '  set',name,'-->',alias
            return
        if self.verbose: print 'Outputs for alias',name,'do not exist:',aliases

    def printStats(self):
        print 'Output       Units    Mean         Min          Max          Stdev'
        print '------------ -------- ------------ ------------ ------------ ------------'
        for output,units in zip(self.outputs,self.units):
            data = getattr(self,output)
            print '{:12s} {:8s} {:12g} {:12g} {:12g} {:12g}'.format(
                        output,
                        units,
                        np.min(data),
                        np.max(data),
                        np.mean(data),
                        np.std(data)
                    )
        print ''

    def plot(self,outputName,*args,**kwargs):
        data = getattr(self,outputName)
        plt.plot(self.t, data, *args, **kwargs)

