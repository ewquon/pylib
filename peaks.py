#!/usr/local/bin/python
import numpy as np

DEBUG       = False
DUMP_PEAKS  = False
TIMING      = True

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
    indices = []
    for i in range(1,N-1):
        if not np.sign(dydx[i]) == np.sign(dydx[i+1]):
            if np.abs(dydx[i]) < np.abs(dydx[i+1]):
                indices.append(i)
            else: indices.append(i+1)

    signchange = np.sign(dydx[1:]) + np.sign(dydx[:N-1])
    indices2 = np.nonzero(signchange==0)[0]
    idxshift = np.abs(dydx[indices2+1]) - np.abs(dydx[indices2]) # < 0 ==> dydx[i+1] closer to 0
    idxshift = (idxshift < 0).astype(int)
    #print np.array(indices)
    #print idxshift
    #print indices2+idxshift
    indices2 += idxshift # correct peak locations
    exact0s = np.nonzero( dydx==0 )[0]
    print indices2.shape, exact0s.shape
    indices2 = set(np.concatenate((indices2,exact0s)))
    print 'num indices',len(indices),len(indices2)
    print 'num set(indices)',len(set(indices)),len(set(indices2))
    diff = set(indices).difference(indices2)
    print 'diff',diff
    assert(len(diff)==0)
    for d in diff:
        print d,dydx[d:d+2]

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
