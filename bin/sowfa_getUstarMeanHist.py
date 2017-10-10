#!/usr/bin/env python
import sys
import numpy as np


makeplots = False
for arg in sys.argv[1:]:
    if arg.lower().strip() == '-plot':
        import matplotlib.pyplot as plt
        makeplots = True
    else:
        logfile = arg

nsteps = 0
startTime = -1
simTimes = []
uStarMean = []
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
except NameError:
    sys.exit('USAGE: '+sys.argv[0]+' log_file')
except IOError:
    sys.exit('Problem reading '+logfile)


np.savetxt('uStarMeanHist.dat',np.vstack((simTimes,uStarMean)).T,fmt='%g',delimiter='\t')


if makeplots:
    plt.figure()
    plt.plot(simTimes,uStarMean)
    plt.xlabel('time step [s]')
    plt.ylabel('uStarMean [m/s]')
    plt.show()
