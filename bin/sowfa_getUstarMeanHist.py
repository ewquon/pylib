#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import numpy as np

histfile = 'uStarMeanHist.dat'

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
uStarMean = []
restarts = []
for logfile in logfiles:
    print('Processing',logfile)
    try:
        with open(logfile,'r') as f:
            for line in f:
                if line.startswith('Create mesh'):
                    startTime = float(line.split()[-1])
                    if startTime > 0:
                        restarts.append(startTime)
                        print('Detected restart from t =',startTime)
                elif line.startswith('Time ='):
                    curTime = float(line.split()[2])
                    simTimes.append(curTime)
                    nsteps += 1
                elif line.startswith('uStarMean ='):
                    us = float(line.split()[2])
                    try:
                        uStarMean[nsteps-1] = us
                    except IndexError:
                        uStarMean.append(us)
    except IOError:
        sys.exit('Problem reading '+logfile)

# prune repeated parts of the time series
simTimes = np.array(simTimes)
uStarMean = np.array(uStarMean)
for rst in restarts:
    rstTime = simTimes[simTimes > rst][0]
    i = np.nonzero(simTimes == rstTime)[0]
    assert(len(i)==2)
    simTimes = np.concatenate((simTimes[:i[0]], simTimes[i[1]:]))
    uStarMean = np.concatenate((uStarMean[:i[0]], uStarMean[i[1]:]))

# if simulation still running, array lengths may differ
Nt,Nu = len(simTimes),len(uStarMean)
if not Nt == Nu:
    N = min(Nt,Nu)
    simTimes = simTimes[:N]
    uStarMean = uStarMean[:N]

print('Writing',histfile)
np.savetxt(histfile, np.vstack((simTimes,uStarMean)).T, fmt='%g', delimiter='\t')


if makeplots:
    plt.figure()
    plt.plot(simTimes,uStarMean)
    plt.xlabel('time step [s]')
    plt.ylabel('uStarMean [m/s]')
    plt.show()
