#!/usr/bin/python
import sys,os
import numpy as np

def pretty_list(strlist,indent=2,sep='\t',width=80):
    """For formatting long lists of strings of arbitrary length
    """
    sep = sep.expandtabs()
    max_item_len = max([len(s) for s in strlist])
    items_per_line = (width - (indent+max_item_len)) / (len(sep)+max_item_len) + 1
    Nlines = len(strlist) / items_per_line
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
    """For post-processing simple xyz-style output from uniform (line) set sample
    """

    sampleExt = 'xy'

    def __init__(self,path='.'):
        """Sets timeNames and sampleNames with the time directory and sample names, respectively.
        Sets t to a numpy array of times
        """
        self.path = path
        self.t = []
        self.timeNames = []
        self.sampleNames = []

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
        self.N = len(self.t)
        if self.N==0: 
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
            assert( fsplit[-1] == self.sampleExt )
            name = '.'.join( fsplit[:-1] )
            self.sampleNames.append( name )

    def __repr__(self):
        s = 'Read times from {:s} :\n{:s}\n'.format(self.path,self.t) \
          + 'Sample names :\n' \
          + pretty_list(sorted(self.sampleNames))
        return s

    def getSample(self,name,field,verbose=True):
        """Returns a sampled field with the specified field name,
        assuming 'x' is identical for all samples

        For a scalar field, the output array has shape (N,NX), where NX
        is the number of spatial samples and N is the length of the time
        array; for a vector field, the output array has shape (N,NX,3).
        """
        found = False
        suffix = '_' + field
        for f in self.sampleNames:
            if f.startswith(name) and f.endswith(suffix):
                found = True
                break
        if not found:
            print 'Sample',name,'with field',field,'not found'
            return

        fname = f + '.' + self.sampleExt

        xfile = self.path + os.sep + f + '.x.npy'
        ufile = self.path + os.sep + f + '.'+field + '.npy'
        try:
            # read from pre-processed numpy data file
            x = np.load( xfile )
            U = np.load( ufile )
            self.NX = U.shape[1]
            print 'Data read from',ufile

        except IOError: # default operation

            # get position of field
            pos = -(len(f) - f.index(suffix))/2
        
            # get number of lines (i.e., number of points in sample) and sampling locations
            testfile = os.path.join(self.path, self.timeNames[0], fname)
            with open(testfile, 'r') as f:
                for i,_ in enumerate(f): pass
            self.NX = i+1

            print 'Found set in',fname,': field',field,'at position',pos,'with',self.NX,'samples'

            x = np.zeros(self.NX)
            with open(testfile, 'r') as f:
                for i,line in enumerate(f):
                    x[i] = float( line.split()[0] )

            # allocate
            if field=='U': 
                isVector = True
                U = np.zeros( (self.N, self.NX, 3) )
            else: 
                isVector = False
                U = np.zeros( (self.N, self.NX) )
            
            # process files in all time dirs
            for it,tdir in enumerate(self.timeNames):
                sys.stdout.write('\r  reading t= %f' % self.t[it])
                with open(os.path.join(self.path, tdir, fname), 'r') as f:
                    if isVector:
                        for i,line in enumerate(f):
                            #Ux[it,i] = float( line.split()[3*pos+0] )
                            #Uy[it,i] = float( line.split()[3*pos+1] )
                            #Uz[it,i] = float( line.split()[3*pos+2] )
                            U[it,i,:] = [ float(val) for val in line.split()[3*pos:][:3] ]
                    else:
                        for i,line in enumerate(f):
                            U[it,i] = float( line.split()[pos] )
            print ''

            print '  saving',ufile
            try:
                np.save( xfile, x )
                np.save( ufile, U )
            except IOError as err:
                print '  warning, unable to write out npy file:',err

        # done reading sample
        if verbose:
            print '  x min/max : [',np.min(x),np.max(x),']'
            if len(U.shape)==3:
                #if isVector:
                print ' ',field+'x min/max : [',np.min(U[:,:,0]),np.max(U[:,:,0]),']'
                print ' ',field+'y min/max : [',np.min(U[:,:,1]),np.max(U[:,:,1]),']'
                print ' ',field+'z min/max : [',np.min(U[:,:,2]),np.max(U[:,:,2]),']'
            else:
                print ' ',field+' min/max : [',np.min(U[:,:]),np.max(U[:,:]),']'

        return x,U


#===============================================================================
#===============================================================================
#===============================================================================
# test run
if __name__=='__main__':

    line0 = uniform(path='linesTransverse')
    print line0

    x,U = line0.getSample('line_08km','U')


