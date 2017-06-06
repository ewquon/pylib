#!/usr/bin/env python
#
# Extract the average wall heat flux from the log file
#
import sys

import refs

if len(sys.argv) <= 1:
    sys.exit('specify log file(s) to process')

ref = refs.read({'rho':1.2, 'Cp':1000.0})

t = []
qwMean = []
for fname in sys.argv[1:]:
    print 'Processing',fname
    with open(fname,'r') as f:
        for line in f:
            if line.startswith('Time ='):
                t.append(float(line.split()[2]))
            elif 'qwMean' in line:
                findstr = 'qwMean ='
                val = line[line.find(findstr)+len(findstr):].split()[0]
                qwMean.append(float(val))

numPredictorCorrectors = len(qwMean)/len(t)
print 'Detected number of predictor and correctors =',numPredictorCorrectors

# pull out the value from the last corrector step
qwMean = qwMean[numPredictorCorrectors-1:-1:numPredictorCorrectors]
qwMean_Wpm2 = qwMean[-1] * ref.rho * ref.Cp # convert to W/m^2
print 'qwMean(t={:.1f} s) = {:.4g} K-m/s = {:.4g} W/m^2'.format(t[-1],qwMean[-1],qwMean_Wpm2)

# save extracted data
import numpy as np
np.savetxt('postProcessing/qwMeanHist.dat',np.array((t,qwMean)).T,fmt='%.4g')

# plot history
import matplotlib.pyplot as plt
plt.plot(t,qwMean)
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$q_w$ [W/m^2]')
plt.suptitle('Average Heat Flux')

plt.savefig('postProcessing/qwMeanHist.png')
plt.show()
