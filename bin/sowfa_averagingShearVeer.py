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
ref = refs.read({'zref':90.0,'Uref':8.0})

data = read( *sys.argv[1:], varList=['U_mean','V_mean'] )

Hhub = ref.zref
heights = [Hhub/4, Hhub/2, Hhub]

Uinf = ref.Uref
alpha = data.shear(heights=heights,Uref=Uinf,zref=Hhub)
veer = data.veer(hmax=Hhub)
print 'shear coefficient, alpha=',alpha
print 'veer=',veer,'deg'

fig,ax = plt.subplots(ncols=2)
ax[0].plot(data.Uh,data.hLevelsCell,'b',label='spatial average')
ax[0].plot(data.approxU,data.approxHeights,'ko',label='fit data')
ax[0].plot(data.approxWindProfile,data.hLevelsCell,'k--',label='curve fit')
ax[0].set_ylabel('height [m]')
ax[0].set_xlabel('horizontal wind [m/s]')
ax[1].plot(data.windDir,data.hLevelsCell)
ax[1].set_xlabel('wind direction [deg]')
ax[0].legend(loc='best')
plt.show()
