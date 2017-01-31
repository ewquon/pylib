#!/usr/bin/env python
import sys,os
import numpy as np
import matplotlib.pyplot as plt

sources = ['UX','UY','UZ','T']
plotHeights = [30,60,90,120,250,500]

if len(sys.argv) <= 1:
    sys.exit('specify postProcessing/SourceHistory directory')

prefix = sys.argv[1]
def srcFilePath(src): return os.path.join(prefix,'Source'+src+'History')

for src in sources:
    if not os.path.isfile( srcFilePath(src) ):
        print 'source history files not found, did you mean...'
        for d in os.listdir(prefix):
            dname = os.path.join(prefix,d)
            if os.path.isdir(dname):
                print ' ',dname
        sys.exit()

def readSourceFile(src):
    z = None
    t = []
    dt = []
    S = []
    readStr = ''
    fname = srcFilePath(src)
    with open(fname,'r') as f:
        hline = f.readline().split()
        readStr += ' '.join(hline[0:2])
        z = [ float(val) for val in hline[2:] ]
        readStr += f.readline().strip()
        for line in f:
            vals = [ float(val) for val in line.split() ]
            t.append( vals[0] )
            dt.append( vals[1] )
            Svals = vals[2:]
            if len(Svals)>1: assert(len(Svals)==len(z))
            S.append(Svals)
    print 'Read',readStr,'from',fname
    return np.array(z),t,dt,np.array(S)

fig,ax = plt.subplots(nrows=len(sources)+1,sharex=True)
for isrc,src in enumerate(sources):
    
    z,t,dt,S = readSourceFile(src)

    if isrc==0:
        ax[0].plot(t,dt,'k')
        ax[0].set_ylabel(r'$\Delta t$')

    if len(S[0,:]) > 1: # we actually used computed sources...
        for ip,pltz in enumerate(plotHeights):
            iz = np.argmin(np.abs(z-pltz))
            print 'Plotting source',src,'at',z[iz]
            r = float(ip)/(len(plotHeights)-1)
            ax[isrc+1].plot(t,S[:,iz],color=[1-r,0,0])
    ax[isrc+1].set_ylabel(r'{:s}'.format(src))

ax[-1].set_xlabel(r't')

plt.tight_layout()
plt.savefig('SourceHistories.png')
plt.show()
