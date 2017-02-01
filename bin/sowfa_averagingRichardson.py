#!/usr/bin/env python
#
# Estimate the Richardson number from the averaged T profile
# - should range between -O(0.01) to +O(0.1) for unstable to stable
#
import sys
from SOWFA.postProcessing.averaging import read
import numpy as np
import matplotlib.pyplot as plt

import refs
ref = refs.read({'g':9.81,'zref':90.0,'D':126.0})

data = read( *sys.argv[1:], varList=['T_mean','U_mean','V_mean'] )

z = data.hLevelsCell
T = data.T_mean[-1,:] # calculate for last 
U = np.sqrt( data.U_mean[-1,:]**2 + data.V_mean[-1,:]**2 )
spacings = z[1:3] - z[0:2]
assert( spacings[0] == spacings[1] )
dz = spacings[0]

rotorTop = ref.zref + ref.D/2
idxs = z<=rotorTop
Tmean = np.mean( T[idxs] )

centralFormula = np.array([-1,0,1])/(2*dz) # second-order accurate
dTdz1 = centralFormula.dot(T[:3])
dUdz1 = centralFormula.dot(U[:3])

oneSidedFormula = np.array([-3,4,-1])/(2*dz) # second-order accurate
dTdz0 = oneSidedFormula.dot(T[:3])
dUdz0 = oneSidedFormula.dot(U[:3])

dTdz = (dTdz1-dTdz0)/dz * (-z[0]) + dTdz0
dUdz = (dUdz1-dUdz0)/dz * (-z[0]) + dUdz0

print 'mean T :',Tmean

print 'dT/dz at z=',z[1],':',dTdz1,' (finite difference)'
print 'dT/dz at z=',z[0],':',dTdz0,' (finite difference)'
print 'dT/dz at z=0:',dTdz,' (extrapolated)'

print 'dU/dz at z=',z[1],':',dUdz1,' (finite difference)'
print 'dU/dz at z=',z[0],':',dUdz0,' (finite difference)'
print 'dU/dz at z=0:',dUdz,' (extrapolated)'

Ri = ref.g/Tmean * dTdz / dUdz**2
print 'Ri :',Ri

fig,ax = plt.subplots(ncols=2)
ax[0].plot(U[idxs],z[idxs],'b-',label='averaging output')
ax[0].plot(U[:3],z[:3],'ko',label='finite-diff points')
ax[0].set_xlabel('U [m/s]')
ax[0].set_ylabel('z [m]')

ax[1].plot(T[idxs],z[idxs],'b-',label='averaging output')
ax[1].plot(T[:3],z[:3],'ko',label='finite-diff points')
ax[1].set_xlabel('T [K]')
ax[1].set_ylabel('z [m]')
ax[1].legend(loc='best')
plt.show()
