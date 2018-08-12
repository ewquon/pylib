#!/usr/bin/env python
#
# This naively distributes sampling planes uniformly through the domain
# based on a provided estimate of the longitudinal and lateral length
# scales (specified separation distance). The grid of sampling planes is
# rotated based on the specified mean wind, and planes that lie outside
# the domain are removed.
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ds = 10.0
Lx = 4320.0         # SOWFA domain length
Ly = 4320.0         # SOWFA domain width
zhub = 90.0         # reference height
L_turbsim = 145.0   # TurbSim grid height/width
Umean,Vmean = 2.2651, 7.67294  # mean velocity components
Sx = 1000.0         # longitudinal spacing between sampling planes
Sy = 500.0          # lateral spacing between sampling planes

wdir = np.arctan2(-Umean,-Vmean)
if wdir < 0: wdir += 2*np.pi
sys.stderr.write('mean wind dir: {:f} deg\n'.format(180./np.pi*wdir))

mag = np.sqrt(Umean**2 + Vmean**2)
nvec1 = (Umean/mag, Vmean/mag, 0)
nvec2 = (-nvec1[1], nvec1[0], 0)
e1 = str(nvec1).replace(',','')
e2 = str(nvec2).replace(',','')

header = """
      inflowPlanes
      {
          type                  sets;
          functionObjectLibs    ("libsampling.so" "libuserfileFormats.so");
          enabled               true;
          interpolationScheme   cellPoint;
          outputControl         timeStep;
          outputInterval        1;
          setFormat             vtkStructured;
          fields
          (
              U
          );
          sets
          ("""

footer = """
          );
      }
"""

inflow_entry = """
              inflowPlane{id:s}
              {{
                  type                array;
                  axis                xyz;
                  origin              ( {x0:g} {y0:g} {z0:g} ); // sampling starts at origin+delta, not origin
                  coordinateRotation
                  {{
                        type          axesRotation;
                        e1            """+e1+""";
                        e2            """+e2+""";
                  }}
                  // TODO: change sampling extent
                  pointsDensity       ( 1 {Ny:d} {Nz:d} );
                  spanBox             (20 {Ly:g} {Lz:g} );
              }}"""
#print(inflow_entry)

# calculate turbsim rotor plane extents
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

# calculate spacing between inflow planes
Sx = np.ceil(Sx/ds) * ds
Sy = np.ceil(Sy/ds) * ds

fig,ax = plt.subplots()
ax.add_patch(
    Rectangle((0,0), Lx, Ly, color='k', lw=3, fill=False)
)

# uniformly layout the plane centers
xr = np.arange(Sx/2,Lx-Sx/2,Sx)
yr = np.arange(Sy/2,Ly-Sy/2,Sy)
xc,yc = np.meshgrid(xr,yr,indexing='ij')
xc = xc.ravel()
yc = yc.ravel()
ax.plot(xc,yc,'kx')
coords = np.stack((xc - Lx/2, yc - Ly/2))
ang = 1.5*np.pi - wdir
Rmat = np.array([[np.cos(ang),-np.sin(ang)],
                 [np.sin(ang), np.cos(ang)]])
coords = np.matmul(Rmat,coords)
coords[0,:] += Lx/2
coords[1,:] += Ly/2
ax.plot(coords[0,:],coords[1,:],'ko') # rotated

# setup rotated samplign planes (and remove those that lie outside the domain)
def in_domain(xp,yp):
    if (np.min(xp) < ds/2) or (np.min(yp) < ds/2) \
            or (np.max(xp) > Lx-ds/2) or (np.max(yp) > Ly-ds/2):
        return False
    else:
        return True
Rmat = np.array([nvec1[:2],nvec2[:2]]).T

selected = []
extents = np.stack([[0,0],yextent])
sys.stderr.write(str(extents)+'\n')
extents = np.matmul(Rmat, extents)
sys.stderr.write(str(extents)+'\n')
for i in range(len(xc)):
    xp = extents[0,:] + coords[0,i]
    yp = extents[1,:] + coords[1,i]
    ax.plot(xp,yp,color='0.8')
    if in_domain(xp,yp):
        sys.stderr.write('plane {:d} extents: {} {}\n'.format(len(selected),xp,yp))
        selected.append(i)

# need workaround because OpenFOAM array sampling only rotates about (0,0,0)
#coordsR = coords[:,selected]
#coordsR[1,:] += yextent[0]
coordsR = coords[:,selected]
coordsR[0,:] += extents[0,0]
coordsR[1,:] += extents[1,0]
for i,plane in enumerate(coordsR.T):
    # add text labels before rotating
    ax.text(plane[0],plane[1],str(i),color='r')
coordsR = np.matmul(Rmat.T, coordsR)

# now create list of inflow planes
inflowplanes = np.stack((coordsR[0,:],
                         coordsR[1,:],
                         zextent[0]*np.ones(len(selected)))) # shape==(3,N)
inflowplanes = inflowplanes.T # shape==(N,3) for enumeration later

ax.plot(coords[0,selected], coords[1,selected], 'bo') # highlight selected

# dump sampling file to stdout
print(header)
for i,plane in enumerate(inflowplanes):
    print(
        inflow_entry.format(
            id='_{:02d}'.format(i),
            x0=plane[0]-ds,
            y0=plane[1]-ds,
            z0=plane[2]-ds,
            Ny=Ny,
            Nz=Nz,
            Ly=(Ny+1)*ds,
            Lz=(Nz+1)*ds,
        )
    )
print(footer)

plt.axis('equal')
plt.show()
