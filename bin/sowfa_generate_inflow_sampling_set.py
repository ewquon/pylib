#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ds = 10.0
zhub = 90.0         # reference height
Lx = 6000.0         # SOWFA domain length
Ly = 3000.0         # SOWFA domain width
L_turbsim = 145.0   # TurbSim grid height/width
L_iso = 135.0       # integral length scale

xoutlet = Lx - 9.5*ds # location of first row
Nrows = 3 # number of rows of sampling planes to generate

inflow_entry = """
              inflowPlane{id:s}
              {{
                  type                array;
                  axis                xyz;
                  origin              ( {x0:g} {y0:g} {z0:g} ); // sampling starts at origin+delta, not origin
                  coordinateRotation
                  {{
                        type          axesRotation;
                        e1            (1 0 0);
                        e2            (0 1 0);
                  }}
                  // TODO: change sampling extent
                  pointsDensity       ( 1 {Ny:d} {Nz:d} );
                  spanBox             (20 {Ly:g} {Lz:g} );
              }}"""

# calculate turbsim extents
yextent = np.array([-L_turbsim/2, L_turbsim/2])
zextent = zhub + yextent
sys.stderr.write('Turbsim y extent: '+str(yextent)+'\n')
sys.stderr.write('Turbsim z extent: '+str(zextent)+'\n')

# round up/down to nearest cell centers
yextent[0] = np.floor(yextent[0]/(ds/2.0)) * ds/2.0
yextent[1] =  np.ceil(yextent[1]/(ds/2.0)) * ds/2.0
zextent[0] = np.floor(zextent[0]/(ds/2.0)) * ds/2.0
zextent[1] =  np.ceil(zextent[1]/(ds/2.0)) * ds/2.0
Ny = int((yextent[1] - yextent[0])/ds) + 1
Nz = int((zextent[1] - zextent[0])/ds) + 1
sys.stderr.write('Sampling y extent: '+str(yextent)+' '+str(Ny)+' pts\n')
sys.stderr.write('Sampling z extent: '+str(zextent)+' '+str(Nz)+' pts\n')

# figure out how many inflow sampling planes to put in a line
Lref = np.ceil(L_iso/ds) * ds
Lref10 = np.ceil(10.0*L_iso/ds) * ds
Nline = int(Ly / (L_turbsim + Lref))
sys.stderr.write('Sampling planes per row: '+str(Nline)+'\n')

# layout planes between Lref0 and Ly-Lref0
x = xoutlet
isim = 0
fig,ax = plt.subplots()
ax.add_patch(
    Rectangle((0,0), Lx, Ly, color='k', lw=3, fill=False)
)
for irow in range(Nrows):
    sys.stderr.write('row'+str(irow)+' : x ='+str(x)+'\n')
    ycenter = np.mod(irow,2)*Lref \
            + (L_turbsim + Lref)*(np.arange(Nline) + 0.5)
    for i,yc in enumerate(ycenter):
        yrange = yc + yextent
        sys.stderr.write('plane'+str(isim)+': '+str(yrange)+'\n')
        print(
            inflow_entry.format(id='_{:02d}'.format(isim),
                                x0=x-ds,
                                y0=yrange[0]-ds,
                                z0=zextent[0]-ds,
                                Ny=Ny,
                                Nz=Nz,
                                Ly=(Ny+1)*ds,
                                Lz=(Nz+1)*ds,
                                )
        )
        ax.plot((x,x), yrange)
        isim += 1
    x -= Lref10

plt.axis('equal')
plt.show()
