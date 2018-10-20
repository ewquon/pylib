#!/usr/bin/env python
from __future__ import print_function
import glob
import sys
import numpy as np

histfile = 'moveDynamicMeshHist.dat'

makeplots = False
logfiles = []
for arg in sys.argv[1:]:
    if arg.lower().strip() == '-plot':
        import matplotlib.pyplot as plt
        makeplots = True
    else:
        logfiles.append(arg)
if len(logfiles) == 0:
    logfiles = glob.glob('log.*Solver')
if len(logfiles) == 0:
    sys.exit('USAGE: '+sys.argv[0]+' logfile(s)')

logfiles.sort()

nsteps = 0
startTime = -1
simTimes = []
restarts = []

# quantities to scrape
maxAR = []
maxNonOrtho = []
avgNonOrtho = []
maxSkewness = []

for logfile in logfiles:
    print('Processing',logfile)
    try:
        with open(logfile,'r') as f:
            for line in f:
                line = line.lstrip()
                if line.startswith('Create mesh'):
                    startTime = float(line.split()[-1])
                    if startTime > 0:
                        restarts.append(startTime)
                        print('Detected restart from t =',startTime)
                elif line.startswith('Time ='):
                    curTime = float(line.split()[2])
                    simTimes.append(curTime)
                    nsteps += 1
                elif line.startswith('Max aspect ratio'):
                    maxAR.append(float(line.split()[-2]))
                elif line.startswith('Mesh non-orthogonality'):
                    line = line.split()
                    maxNonOrtho.append(float(line[-3]))
                    avgNonOrtho.append(float(line[-1]))
                elif line.startswith('Max skewness'):
                    maxSkewness.append(float(line.split()[-2]))
    except IOError:
        sys.exit('Problem reading '+logfile)

# prune repeated parts of the time series
simTimes = np.array(simTimes)
maxAR = np.array(maxAR)
maxNonOrtho = np.array(maxNonOrtho)
avgNonOrtho = np.array(avgNonOrtho)
maxSkewness = np.array(maxSkewness)
for rst in restarts:
    rstTime = simTimes[simTimes > rst][0]
    i = np.nonzero(simTimes == rstTime)[0]
    if len(i) == 2:
        print('Pruning time series at t= {:f} s'.format(rstTime))
        simTimes = np.concatenate((simTimes[:i[0]], simTimes[i[1]:]))
        maxAR = np.concatenate((maxAR[:i[0]], maxAR[i[1]:]))
        maxNonOrtho = np.concatenate((maxNonOrtho[:i[0]], maxNonOrtho[i[1]:]))
        avgNonOrtho = np.concatenate((avgNonOrtho[:i[0]], avgNonOrtho[i[1]:]))
        maxSkewness = np.concatenate((maxSkewness[:i[0]], maxSkewness[i[1]:]))

# if simulation still running, array lengths may differ
Nt,N = len(simTimes),len(maxSkewness)
if not Nt == N:
    N = min(Nt,N)
    simTimes = simTimes[:N]
    maxAR = maxAR[:N]
    maxNonOrtho = maxNonOrtho[:N]
    avgNonOrtho = avgNonOrtho[:N]
    maxSkewness = maxSkewness[:N]

print('Writing',histfile)
np.savetxt(histfile,
        np.vstack((simTimes,maxAR,maxNonOrtho,avgNonOrtho,maxSkewness)).T,
        fmt='%g', delimiter=',')


if makeplots:
    fig,ax = plt.subplots(nrows=4,sharex=True)
    ax[0].plot(simTimes,maxAR)
    ax[1].plot(simTimes,maxNonOrtho)
    ax[2].plot(simTimes,avgNonOrtho)
    ax[3].plot(simTimes,maxSkewness)
    ax[0].set_title('max aspect ratio')
    ax[1].set_title('max mesh non-orthogonality')
    ax[2].set_title('average mesh non-orthogonality')
    ax[3].set_title('max skewness')
    ax[-1].set_xlabel('time step [s]')
    plt.show()
