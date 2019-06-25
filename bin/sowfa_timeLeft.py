#!/usr/bin/env python
import sys
import os
import glob
import myutils
import time
import numpy as np

log_est_total_time = 'estimated_total_time'

dpath = '.'

makeplots = False
if '-plot' in sys.argv:
    makeplots = True
    import matplotlib.pyplot as plt
    sys.argv.remove('-plot')

writehist = False
if '-write' in sys.argv:
    writehist = True
    sys.argv.remove('-write')

if len(sys.argv) > 1:
    fpath = sys.argv[1]
    dpath = os.path.split(fpath)[0]

with open(os.path.join(dpath,'system','controlDict'),'r') as f:
    for line in f:
        if line.strip().startswith('application'):
            app = line.split()[1].split(';')[0]
        elif line.strip().startswith('deltaT'): 
            dt = float(line.split()[1][:-1])
        elif line.strip().startswith('endTime'): 
            endTime = float(line.split()[1][:-1])

if len(sys.argv) == 1:
    import glob
    logfiles = glob.glob('log.*'+app)
    logfiles.sort()
    fpath = logfiles[-1]

nsteps = 0
startTime = -1
simTimes = []
clockTimes = []
CourantMax = []
try:
    with open(fpath,'r') as f:
        for line in f:
            if line.startswith('Create mesh'):
                startTime = float(line.split()[-1])
                if startTime > 0:
                    print('Detected restart from t =',startTime)
            elif line.startswith('Time ='):
                curTime = float(line.split()[2])
                simTimes.append(curTime)
                nsteps += 1
            elif line.startswith('ExecutionTime ='):
                elapsedTime = float(line.split()[6])
                clockTimes.append(elapsedTime)
            elif line.startswith('Courant Number'):
                cflmax = float(line.split()[-1])
                CourantMax.append(cflmax)
except NameError:
    sys.exit('USAGE: '+sys.argv[0]+' log_file')
except IOError:
    sys.exit('Problem reading '+fpath)

completed = (curTime-startTime) / (endTime-startTime)

print('Simulation is at currently at t = {:.1f} and will end at {:.1f} ({:.1f}% complete)'.format(
      curTime, endTime, 100*completed))
print(app,'has been running for',nsteps,'steps')
print('Elapsed time:',myutils.smartTime(elapsedTime))

timestep_size = np.diff(np.array(simTimes))
ctime_per_step = np.diff(np.array(clockTimes))
print('Average/min/max time per step', \
    np.mean(ctime_per_step), np.min(ctime_per_step), np.max(ctime_per_step))

totalTime = elapsedTime / completed
print('ESTIMATED TOTAL TIME:',myutils.smartTime(totalTime))
if log_est_total_time:
    try:
        with open(os.path.join(dpath,log_est_total_time),'a') as f:
            f.write(time.strftime('%x %X\tsimulated t={:.1f}s, N={:d}\tESTIMATED TOTAL TIME: {:s}\n'.format(curTime,nsteps,myutils.smartTime(totalTime))))
    except IOError:
        pass

remainingTime = totalTime - elapsedTime
print('Remaining time:',myutils.smartTime(remainingTime))

timeStepsLeft = (endTime-curTime) / timestep_size[-1]
remainingTimeOpt = timeStepsLeft * ctime_per_step[-1]
print('Remaining time (based on clocktime for last timestep):',myutils.smartTime(remainingTimeOpt))
print(' ',timeStepsLeft,'time steps left')
print('  ~',ctime_per_step[-1],'clock time per step')

print('Average/min/max maximum Courant number per step', \
    np.mean(CourantMax), np.min(CourantMax), np.max(CourantMax))

print(time.strftime('Current date/time is %x %X'))

if writehist:
    N = min(len(simTimes),len(ctime_per_step))
    with open(fpath+'.timing','w') as f:
        f.write(time.strftime('# Current system time: %x %X\n'))
        f.write('# step  simulation_time  wallclock_time_per_step\n')
        for i in range(N):
            f.write('{:d} {:f} {:f}\n'.format(i+1,simTimes[i],ctime_per_step[i]))

if makeplots:
    plt.figure()
    plt.plot(ctime_per_step)
    plt.xlabel('iteration')
    plt.ylabel('Clock time / step [s]')

    if (np.max(timestep_size)-np.min(timestep_size)) > 1e-6:
        plt.figure()
        plt.plot(timestep_size)
        plt.xlabel('iteration')
        plt.ylabel('simulated timestep [s]')
    else:
        print('Skipping timestep plot for constant dt=',np.mean(timestep_size))

    plt.show()
