#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

if len(sys.argv)==1:
    sys.exit('USAGE: '+sys.argv[0]+' csvFile [tStart]')

tstart = 0.
if len(sys.argv)>2:
    tstart = int(sys.argv[2])

data = np.loadtxt(sys.argv[1],delimiter=',',skiprows=1)
t = data[:,0]
TI = data[:,1] # convert to percent

fig,ax = plt.subplots(figsize=(4,3))
ax.plot(t,100*TI)
ax.set_xlabel('time [s]')
ax.set_ylabel('TI [%]')

trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
ax.fill_between(t, 0, 1, where=t<tstart, edgecolor='', facecolor='grey', alpha=0.5, transform=trans)

TIfinal = TI[t >= tstart]
print 'TI Statistics for t > {:f} ({:d} data points)'.format(tstart,len(TIfinal))
print 'approx stats window',t[-1]-tstart,'s'
print 'mean  :',np.mean(TIfinal)
print 'stdev :',np.std(TIfinal)

plt.tight_layout()
plt.show()

outname = sys.argv[1].split('.')
outname = '.'.join(outname[:-1]) +'_'+str(tstart) + '.png'
fig.savefig(outname)
print 'Wrote out',outname
