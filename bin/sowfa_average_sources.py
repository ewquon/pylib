#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
T = 1e5

sources = dict()
def process(s):
    s = s.replace('(','')
    s = s.replace(')','')
    s = s.replace(';','')
    s = s.split()
    name = s.pop(0)
    sources[name] = np.array(s, dtype=float)
with open('sources','r') as f:
    curstr = ''
    for line in f:
        curstr += line
        if line.strip().endswith(';'):
            process(curstr)
            curstr = ''
sys.stderr.write('Read: '+str(sources.keys())+'\n')

#
# reshape sources arrays into time-height arrays
#
z = sources['sourceHeightsMomentum']
z_tmp = sources['sourceHeightsTemperature']
assert(np.all(z== z_tmp))
NZ = len(z)
NT = len(z_tmp)
xmom = sources['sourceTableMomentumX']
N = len(xmom) / (NZ+1)
xmom = xmom.reshape((N,NZ+1))
t = xmom[:,0]
xmom = xmom[:,1:]
ymom = sources['sourceTableMomentumY'].reshape((N,NZ+1))[:,1:]
zmom = sources['sourceTableMomentumZ'].reshape((N,NZ+1))[:,1:]
temp = sources['sourceTableTemperature'].reshape((N,NZ+1))[:,1:]

#
# average starting from specified time
#
if len(sys.argv) > 1:
    t0 = float(sys.argv[1])
else:
    t0 = 0.0
Navg = len(t[t>=t0])
sourceU = np.mean(xmom[t>=t0,:], axis=0)
sourceV = np.mean(ymom[t>=t0,:], axis=0)
sourceW = np.mean(zmom[t>=t0,:], axis=0)
sourceT = np.mean(temp[t>=t0,:], axis=0)
sys.stderr.write('Averaging from t={:.1f} to {:.1f} over {:d} iterations\n'.format(t0,t[-1],Navg))

fig,ax = plt.subplots(ncols=4,figsize=(6,6))
ist = np.nonzero(t>=t0)[0][0]
ax[0].plot(sourceU,z)
ax[0].plot(np.min(xmom[t>=t0],axis=0),z,'k--')
ax[0].plot(np.max(xmom[t>=t0],axis=0),z,'k--')
ax[0].set_xlabel('sourceU')
ax[1].plot(sourceV,z)
ax[1].plot(np.min(ymom[t>=t0],axis=0),z,'k--')
ax[1].plot(np.max(ymom[t>=t0],axis=0),z,'k--')
ax[1].set_xlabel('sourceV')
ax[2].plot(sourceW,z)
ax[2].plot(np.min(zmom[t>=t0],axis=0),z,'k--')
ax[2].plot(np.max(zmom[t>=t0],axis=0),z,'k--')
ax[2].set_xlabel('sourceW')
ax[3].plot(sourceT,z)
ax[3].plot(np.min(temp[t>=t0],axis=0),z,'k--')
ax[3].plot(np.max(temp[t>=t0],axis=0),z,'k--')
ax[3].set_xlabel('sourceT')

#
# output new sources file
#
print('sourceHeightsMomentum')
print('(')
for zi in z:
    print('    {:f}'.format(zi))
print(');')
print('\nsourceTableMomentumX')
print('(')
srcstr = ' '.join([str(val) for val in sourceU])
print('    ({:f} {:s})'.format(0,srcstr))
print('    ({:f} {:s})'.format(T,srcstr))
print(');')
print('\nsourceTableMomentumY')
print('(')
srcstr = ' '.join([str(val) for val in sourceV])
print('    ({:f} {:s})'.format(0,srcstr))
print('    ({:f} {:s})'.format(T,srcstr))
print(');')
print('\nsourceTableMomentumZ')
print('(')
srcstr = ' '.join([str(val) for val in sourceW])
print('    ({:f} {:s})'.format(0,srcstr))
print('    ({:f} {:s})'.format(T,srcstr))
print(');')
print('\nsourceHeightsTemperature')
print('(')
for zi in z:
    print('    {:f}'.format(zi))
print(');')
print('\nsourceTableTemperature')
print('(')
srcstr = ' '.join([str(val) for val in sourceT])
print('    ({:f} {:s})'.format(0,srcstr))
print('    ({:f} {:s})'.format(T,srcstr))
print(');')
fig.savefig('source_mean.png')

#
# check convergence
#
#xmom_converge = np.max(np.abs(np.diff(xmom, axis=0)), axis=1)
#ymom_converge = np.max(np.abs(np.diff(ymom, axis=0)), axis=1)
#zmom_converge = np.max(np.abs(np.diff(zmom, axis=0)), axis=1)
#temp_converge = np.max(np.abs(np.diff(temp, axis=0)), axis=1)
xmom_converge = np.std(np.diff(xmom, axis=0), axis=1)
ymom_converge = np.std(np.diff(ymom, axis=0), axis=1)
zmom_converge = np.std(np.diff(zmom, axis=0), axis=1)
temp_converge = np.std(np.diff(temp, axis=0), axis=1)

fig,ax = plt.subplots()
ax.semilogy(t[1:], xmom_converge, label='U')
ax.semilogy(t[1:], ymom_converge, label='V')
ax.semilogy(t[1:], zmom_converge, label='W')
ax.semilogy(t[1:], temp_converge, label='T')
ax.set_xlabel('Time [s]')
ax.set_ylabel('RMS change in source terms')
ax.legend(loc='best')
fig.savefig('source_rms_hist.png')

#T,Z = np.meshgrid(t,z,indexing='ij')
#print(T.shape,Z.shape,xmom.shape)
#plt.contourf(T[::100,:],Z[::100,:],xmom[::100,:])

plt.show()
