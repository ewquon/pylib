#!/usr/local/bin/python
import numpy as np

DEBUG       = False
DUMP_PEAKS  = False
TIMING      = False

if DEBUG: import matplotlib.pyplot as plt
if TIMING: import time

def find_peaks(data,x=[],Nsmoo=1):

    N = len(data)
    i = np.arange(N)
    if len(x)==0: 
        x = i # assume evenly spaced data
    else:
        assert(len(x)==N)
    dx = np.diff(x)

    if TIMING: t0 = time.time()
    ist = i   - Nsmoo
    ind = i+1 + Nsmoo
    ist = np.maximum(0,ist)
    ind = np.minimum(N,ind)
    width = ind - ist
    if TIMING: print ' - walltime [s] to calculate indices:',time.time()-t0

    if TIMING: t0 = time.time()
    ysmoo = data.copy()
    buff = np.zeros((N))
    for ismoo in range(1,Nsmoo+1):
        # add LHS
        buff[:ismoo] = 0.
        buff[ismoo:] = data[:N-ismoo]
        if DEBUG: print -ismoo,buff
        ysmoo += buff
        # add RHS
        buff[:N-ismoo] = data[ismoo:]
        buff[N-ismoo:] = 0.
        ysmoo += buff
        if DEBUG: print ismoo,buff
    ysmoo /= width
    if TIMING: print ' - walltime [s] to smooth data:',time.time()-t0

    if TIMING: t0 = time.time()
    twodx = dx[1:] + dx[:-1]
    dydx = np.zeros((N))
    dydx[1:-1] = (ysmoo[2:] - ysmoo[:-2]) / twodx
    if TIMING: print ' - walltime [s] to calculate slope:',time.time()-t0

    if TIMING: t0 = time.time()
    signs = np.sign(dydx)
    signchange = signs[1:] + signs[:N-1]
    idxchange = np.nonzero(signchange==0)[0]
    idxshift = np.abs(dydx[idxchange+1]) - np.abs(dydx[idxchange]) # < 0 ==> dydx[i+1] closer to 0
    idxchange += (idxshift < 0).astype(int) # correct peak locations
    exact0s = np.nonzero( dydx==0 )[0]
    indices = np.unique( np.concatenate((idxchange,exact0s)) )
    if TIMING: print ' - walltime [s] to find peak indices:',time.time()-t0

    peaks = data[indices]

    if DUMP_PEAKS:
        with open('peaks.tmp','w') as f:
            for xi,yi in zip(x[indices],peaks):
                f.write('%f %f\n'%(xi,yi))
        print 'wrote peaks.tmp'

    if DEBUG:
        print len(indices),'min/maxima found',indices
        plt.plot(x,data,label='input data')
        plt.plot(x,ysmoo,'k--',label='smoothed input')
        plt.plot(x[indices],peaks,'o',label='detected peaks')
        #plt.xlim((72,77))
        #plt.ylim((200,700))
        plt.legend(loc='best')
        plt.show()

    return peaks, indices
