#!/usr/bin/env python
#
# Estimate the Richardson number from the averaged T profile at the final time step
# - should range between -O(0.01) to +O(0.1) for unstable to stable
#
import sys
from SOWFA.postProcessing.averaging import read
import numpy as np
import matplotlib.pyplot as plt

import refs
ref = refs.read({'g':9.81,'zref':90.0,'D':126.0})

data = read( *sys.argv[1:], varList=['T_mean','U_mean','V_mean'] )

Ri = data.calcRichardsonNumber(g=ref.g,zref=ref.zref,D=ref.D)
print 'Ri :',Ri

#fig,ax = plt.subplots(ncols=2)
#ax[0].plot(U[idxs],z[idxs],'b-',label='averaging output')
#ax[0].plot(U[:3],z[:3],'ko',label='finite-diff points')
#ax[0].set_xlabel('U [m/s]')
#ax[0].set_ylabel('z [m]')
#
#ax[1].plot(T[idxs],z[idxs],'b-',label='averaging output')
#ax[1].plot(T[:3],z[:3],'ko',label='finite-diff points')
#ax[1].set_xlabel('T [K]')
#ax[1].set_ylabel('z [m]')
#ax[1].legend(loc='best')
#plt.show()
