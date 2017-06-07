#!/usr/bin/python
import sys
import myutils
import time
import numpy as np

with open('system/controlDict','r') as f:
    for line in f:
	if line.strip().startswith('application'):
	    app = line.split()[1].split(';')[0]
        elif line.strip().startswith('deltaT'): 
            dt = float(line.split()[1][:-1])
        elif line.strip().startswith('endTime'): 
            endTime = float(line.split()[1][:-1])

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
clockTimes = []
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
            elif line.startswith('ExecutionTime ='):
                elapsedTime = float(line.split()[6])
                clockTimes.append(elapsedTime)
except IOError:
    sys.exit('USAGE: '+sys.argv[0]+' log_file')

completed = (curTime-startTime) / (endTime-startTime)

print 'Simulation is at currently at t = %.1f and will end at %.1f (%.1f%% complete)' \
    % (curTime, endTime, 100*completed)
print app,'has been running for',nsteps,'steps'
print 'Elapsed time:',myutils.smartTime(elapsedTime)

timestep_size = np.diff(np.array(simTimes))
ctime_per_step = np.diff(np.array(clockTimes))
print 'Average/min/max time per step', \
    np.mean(ctime_per_step), np.min(ctime_per_step), np.max(ctime_per_step)

totalTime = elapsedTime / completed
print 'ESTIMATED TOTAL TIME:',myutils.smartTime(totalTime)

remainingTime = totalTime - elapsedTime
print 'Remaining time:',myutils.smartTime(remainingTime)
print time.strftime('Current date/time is %x %X')

if makeplots:
    plt.figure()
    plt.plot(ctime_per_step)
    plt.xlabel('iteration')
    plt.ylabel('Clock time / step [s]')

    plt.figure()
    plt.plot(timestep_size)
    plt.xlabel('iteration')
    plt.ylabel('simulated timestep [s]')

    plt.show()