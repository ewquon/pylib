#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import numpy as np

histfile = 'BoundaryFluxHist.dat'

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
hist = []
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
                elif line.startswith('Total Boundary Flux'):
                    line = line.split()
                    fluxerr = float(line[-1])
                    try:
                        hist[nsteps-1] = fluxerr
                    except IndexError:
                        hist.append(fluxerr)
    except IOError:
        sys.exit('Problem reading '+logfile)

if len(simTimes) > len(hist):
    simTimes = simTimes[:len(hist)]

if makeplots:
    plt.figure()
    plt.semilogy(simTimes,np.abs(hist))
    plt.xlabel('time step [s]')
    plt.ylabel('Boundary Flux Error')
    plt.show()

print('Writing',histfile)
np.savetxt(histfile, np.vstack((simTimes,hist)).T, fmt='%g', delimiter='\t')

