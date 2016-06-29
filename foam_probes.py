#!/usr/local/bin/python
import numpy as np

class probedata:
    """ Class to process OpenFOAM probe data
    """

    def __init__(self,f):
        self._read_probe_positions(f)
        self._read_data(f)

    def _read_probe_positions(self,f):
        self.pos = []
        line = f.readline()
        while '(' in line and ')' in line:
            line = line.strip()
            assert(line[0]=='#')
            assert(line[-1]==')')
            iprobe = int(line.split()[2])
            i = line.find('(')
            pt = [ float(val) for val in line[i+1:-1].split() ]
            self.pos.append( np.array(pt) )
            line = f.readline()
        self.N = len(self.pos)
        assert(self.N == iprobe+1)
        print self.N, 'probes read'

    def _read_data(self,f):
        line = f.readline()
        assert(line.split()[1] == 'Time')
        t = []
        U = []
        for i in range(self.N): U.append([])
        for line in f:
            t.append(float(line.split()[0]))
            line = line.replace('(',',').replace(')','')
            for i,vec in enumerate(line.split(',')[1:]):
                U[i].append( np.array([ float(val) for val in vec.split() ]) )
        self.t = np.array(t)
        self.U = np.array(U)
        self.Nt = len(self.t)
        print 'Times read:',self.Nt,self.t
        print 'U size:',self.U.shape

    def calc_Umag(self):
        self.Umag = np.zeros((self.N,self.Nt))
        for iprobe in range(self.N):
            for i in range(self.Nt):
                U = self.U[iprobe,i,:]
                self.Umag[iprobe,i] = np.sqrt(U.dot(U))

