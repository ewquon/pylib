#!/usr/bin/env python
import sys
from SOWFA.postProcessing.averaging import read
import matplotlib.pyplot as plt
import numpy as np
heights = [90.]
tavg_window = 600.
dt = 1.0
SFS = True

data = read( *sys.argv[1:] )
tavg,TIx_hist,TIy_hist,TIz_hist,TIdir_hist,TIxyz_hist,TKE_hist = data.calcTI_hist( heights, tavg_window, dt, SFS )
for ih,z in enumerate(heights):
    plt.plot( tavg, TIdir_hist[:,ih], label='z={:.1f}'.format(z) )
    fname = 'TIhist_z{:.1f}.csv'.format(z)
    np.savetxt( fname, np.vstack((tavg,TIdir_hist[:,ih])).T, delimiter=',', header='Time,TI' )
    print 'wrote',fname

plt.legend(loc='best')
plt.show()

plt.savefig('TIhist.png')

data.calcTI( heights )
for i,h in enumerate(heights):
    print 'z=',h,':'
    print '  TIx/y/z =',data.TIx[i],data.TIy[i],data.TIz[i]
    print '  TIdir   =',data.TIdir[i]
    print '  TIxyz   =',data.TIxyz[i],' ( TKE=',data.TKE[i],')'

