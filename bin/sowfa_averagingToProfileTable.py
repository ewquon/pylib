#!/usr/bin/env python
# run sowfa_averagingAll.py first to create averagingProfiles.csv
import sys
import numpy as np

if len(sys.argv) < 3:
    sys.exit('USAGE: '+sys.argv[0]+' [path/to/averagingProfiles.csv] [constant/initProfileTable]')

# header: z,U,V,W,T,uu,vv,ww,uv,uw,vw,R11,R22,R33,R12,R13,R23
tableData = np.loadtxt(sys.argv[1],delimiter=',',skiprows=1)
z = tableData[:,0]
U = tableData[:,1]
V = tableData[:,2]
W = tableData[:,3]
T = tableData[:,4]

with open(sys.argv[2],'w') as f:
    f.write('profileTable\n')
    f.write('(\n')
    f.write('//  z\tU\tV\tT\n')
    for zi,Ui,Vi,Ti in zip(z,U,V,T):
        f.write('    ({:g}  {:g} {:g} {:g})\n'.format(zi,Ui,Vi,Ti))
    f.write(');')
