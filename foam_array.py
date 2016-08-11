#!/usr/bin/python
import sys,os
import numpy as np

class ensight:
    """ For post-processing array output in ensight format
    """

    vectorFields = ['U']

    def __init__(self,path='.'):
        """ Sets timeNames to the time directory names
        Sets t to a numpy array of times
        Sets sampleNames to the name of the scalar or vector fields
        """
        print 'Initializing ensight solution object'
        self.path = path
        self.meshRead = False

        self.timeNames = []
        self.t = None
        self.sampleNames = []
        self.fieldNames = []
        self.fieldTypes = dict()

        self.Nt = -1
        self.N = -1

        # get list of times from directory names
        try:
            dirs = next(os.walk(path))[1]
        except StopIteration:
            print 'Path',path,'not found'
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
        self.Nt = len(self.t)
        if self.Nt==0: 
            print 'No time directories found in', path
            return

        # sort based on times
        order = np.argsort(self.t)
        self.t = self.t[order]
        self.timeNames = [ timeNames[i] for i in order ]
        dt = np.diff(self.t)
        if np.max(dt)-np.min(dt) > 1e-14:
            print 'Warning: gaps detected in sampling times'
            self.dt = dt
        else:
            self.dt = dt[0]

        # gather sample names
        path0 = os.path.join(path,timeNames[0])
        flist = [ f for f in os.listdir(path0) if os.path.isfile(os.path.join(path0, f)) ]
        for f in flist:
            fsplit = f.split('.')
            if not fsplit[-1] == 'case': continue
            with open(os.path.join(path0,f),'r') as f:
                for line in f:
                    if line.startswith('scalar') or line.startswith('vector'):
                        line = line.split()
                        varname = line[4]
                        assert( not varname in self.fieldNames )
                        self.fieldNames.append(varname)
                        self.fieldTypes[varname] = line[0]
                    elif line.startswith('number of steps'):
                        assert( int(line.split()[-1])==1 )
                        name = '.'.join( fsplit[:-1] )
                        self.sampleNames.append( name )
                        break

    def __repr__(self):
        s = 'Read %d times, t = %s\n' % (self.Nt,self.t) \
          + 'Sample names :'
        for name in sorted( self.sampleNames ):
            s += '\n  %s' % name
        s += '\nSampled fields :'
        for field in self.fieldNames:
            s += ' ' + field
        return s

    def readMesh(self,fname=None):
        if fname is None:
            fname = os.path.join(self.path,self.timeNames[0],self.sampleNames[0]+'.mesh')
        print 'Reading mesh from',fname
        with open(fname,'r') as f:
            for _ in range(5): f.readline() # header
            assert( int(f.readline())==1 ) # single part in this file
            f.readline()
            f.readline()
            self.N = int(f.readline()) # number of points
            remainingLines = f.readlines()
        self.x = np.array([ float(val) for val in remainingLines[        :  self.N] ])
        self.y = np.array([ float(val) for val in remainingLines[  self.N:2*self.N] ])
        self.z = np.array([ float(val) for val in remainingLines[2*self.N:3*self.N] ])

        print '  x min/max : [',np.min(self.x),np.max(self.x),']'
        print '  y min/max : [',np.min(self.y),np.max(self.y),']'
        print '  z min/max : [',np.min(self.z),np.max(self.z),']'

        self.meshRead = True
        return self.x, self.y, self.z

    def getSample(self,name,field,verbose=True,updateMesh=False,dump=True):
        """ Returns a field from sample with specified name
        For a scalar field, the output array has shape (N,NX), where NX is the number of spatial samples and N is the length of the time array;
        for a vector field, the output array has shape (N,NX,3)
        """
        if name in self.sampleNames and field in self.fieldNames:
            pass
        else:
            print 'Sample',name,'with field',field,'not found'
            return

        # read mesh if necessary
        if self.meshRead is False or updateMesh is True:
            fname = os.path.join(self.path,self.timeNames[0],name+'.mesh')
            self.readMesh(fname=fname)

        # allocate
        if self.fieldTypes[field]=='vector':
            isVector = True
            u = np.zeros( (self.Nt, self.N, 3) )
        else: 
            isVector = False
            u = np.zeros( (self.Nt, self.N) )
        
        print 'Reading',self.fieldTypes[field],'field',field,'from sample',name

        # process files in all time dirs
        suffix = '.000.' + field
        for it,tdir in enumerate(self.timeNames):
            sys.stdout.write('\r  processing t= %f' % self.t[it])
            with open(os.path.join(self.path, tdir, name+suffix), 'r') as f: data = f.readlines()[4:]
            if isVector:
                u[it,:,0] = [ float(val) for val in data[        :  self.N] ]
                u[it,:,1] = [ float(val) for val in data[  self.N:2*self.N] ]
                u[it,:,2] = [ float(val) for val in data[2*self.N:3*self.N] ]
            else:
                u[it,:] = [ float(val) for val in data ]
        print ''

        # done reading sample
        if verbose:
            if isVector:
                print ' ',field+'x min/max : [',np.min(u[:,:,0]),np.max(u[:,:,0]),']'
                print ' ',field+'y min/max : [',np.min(u[:,:,1]),np.max(u[:,:,1]),']'
                print ' ',field+'z min/max : [',np.min(u[:,:,2]),np.max(u[:,:,2]),']'
            else:
                print ' ',field+' min/max : [',np.min(u[:,:]),np.max(u[:,:]),']'

        # dump np.array for quick read next time
        if dump:
            dumpname = name + '.npz'
            np.savez( dumpname, U=u )
            print 'Dumped velocity ndarray to', dumpname

        return u


#===============================================================================
#===============================================================================
#===============================================================================
# test run
#if __name__=='__main__':

