#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import numpy as np

histfile = 'CourantHist.dat'

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
CoMean = []
CoMax = []
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
                elif line.startswith('Courant Number'):
                    line = line.split()
                    cmean = float(line[3])
                    cmax = float(line[5])
                    try:
                        CoMean[nsteps-1] = cmean
                    except IndexError:
                        CoMean.append(cmean)
                    try:
                        CoMax[nsteps-1] = cmax
                    except IndexError:
                        CoMax.append(cmax)
    except IOError:
        sys.exit('Problem reading '+logfile)

print('Writing',histfile)
np.savetxt(histfile, np.vstack((simTimes,CoMean,CoMax)).T, fmt='%g', delimiter='\t')


if makeplots:
    plt.figure()
    plt.plot(simTimes,CoMean)
    plt.plot(simTimes,CoMax,'k--')
    plt.xlabel('time step [s]')
    plt.ylabel('Courant Number')
    plt.show()
