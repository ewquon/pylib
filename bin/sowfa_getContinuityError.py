#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import numpy as np

histfile = 'ContinuityErrorHist.dat'

makeplots = False
logfiles = []
for arg in sys.argv[1:]:
    if arg.lower().strip() == '-plot':
        import matplotlib.pyplot as plt
        makeplots = True
    else:
        logfiles.append(arg)
if len(logfiles) == 0:
    sys.exit('USAGE: '+sys.argv[0]+' logfile(s)')

nsteps = 0
startTime = -1
simTimes = []
errMin = []
errMax = []
errMean = []
for logfile in logfiles:
    print('Processing',logfile)
    try:
        with open(logfile,'r') as f:
            for line in f:
                if line.startswith('Create mesh'):
                    startTime = float(line.split()[-1])
                    if startTime > 0: print('Detected restart from t =',startTime)
                elif line.startswith('Time ='):
                    curTime = float(line.split()[2])
                    simTimes.append(curTime)
                    nsteps += 1
                elif line.startswith('Local Flux Continuity Error'):
                    line = line.split()
                    minval = float(line[-6])
                    maxval = float(line[-4])
                    meanval = float(line[-1])
                    try:
                        errMin[nsteps-1] = minval
                    except IndexError:
                        errMin.append(minval)
                    try:
                        errMax[nsteps-1] = maxval
                    except IndexError:
                        errMax.append(maxval)
                    try:
                        errMean[nsteps-1] = meanval
                    except IndexError:
                        errMean.append(meanval)
    except IOError:
        sys.exit('Problem reading '+logfile)

print('Writing',histfile)
np.savetxt(histfile, np.vstack((simTimes,errMin,errMax,errMean)).T, fmt='%g', delimiter='\t')


if makeplots:
    plt.figure()
    plt.semilogy(simTimes,errMin,label='min')
    plt.semilogy(simTimes,errMax,label='max')
    plt.semilogy(simTimes,errMean,'k--',label='weighted mean')
    plt.xlabel('time step [s]')
    plt.ylabel('Continuity Error')
    plt.legend(loc='best')
    plt.show()
