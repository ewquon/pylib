#!/usr/bin/env python
import sys
from SOWFA.postProcessing.averaging import read
heights = [80.,90.]
data = read( *sys.argv[1:] )
data.calcTI( heights )

for i,h in enumerate(heights):
    print 'z=',h,':'
    print '  TIx/y/z =',data.TIx[i],data.TIy[i],data.TIz[i]
    print '  TIdir   =',data.TIdir[i]
    print '  TIxyz   =',data.TIxyz[i],' ( TKE=',data.TKE[i],')'
