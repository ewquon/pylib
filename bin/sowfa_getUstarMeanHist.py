#!/usr/bin/env python
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
for logfile in logfiles:
    print 'Processing',logfile
    try:
        with open(logfile,'r') as f:
            for line in f:
                if line.startswith('Create mesh'):
                    startTime = float(line.split()[-1])
                    if startTime > 0: print 'Detected restart from t =',startTime
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

if os.path.isfile(histfile):
    np.savetxt(open(histfile,'a'),np.vstack((simTimes,uStarMean)).T,fmt='%g',delimiter='\t')
    print 'Appending to',histfile
else:
    np.savetxt(open(histfile,'w'),np.vstack((simTimes,uStarMean)).T,fmt='%g',delimiter='\t')
    print 'Writing new',histfile


if makeplots:
    plt.figure()
    plt.plot(simTimes,uStarMean)
    plt.xlabel('time step [s]')
    plt.ylabel('uStarMean [m/s]')
    plt.show()
