#!/usr/bin/env python
#
# Tool to convert array sampling data from OpenFOAM into a series of structured
# VTK files.
#
# Written by Eliot Quon (eliot.quon@nrel.gov) - 2017-07-07
#
import sys,os
import numpy as np
#import foamdict as fd
#from evtk.hl import gridToVTK
import VTKwriter

verbose = False

if len(sys.argv) <= 5:
    sys.exit('USAGE: '+sys.argv[0]+' samplingDir [samplingDefinitionFile] [prefix Nx Ny Nz]')
samplingDir = sys.argv[1]
prefix = sys.argv[2]
if os.path.isfile(prefix):
    #samplingdef = fd.dictfile(prefix)
    sys.exit('Processing of smapling file definition not implemented yet')
else:
    Nx = int(sys.argv[3])
    Ny = int(sys.argv[4])
    Nz = int(sys.argv[5])

N = Nx*Ny*Nz

# get sorted list of time directories and times
listing = os.listdir(samplingDir)
times = []
timeDirs = []
for dirname in listing:
    path = os.path.join(samplingDir,dirname)
    if os.path.isdir(path):
        timeDirName = os.path.split(path)[1]
        try:
            times.append(float(timeDirName))  
        except ValueError:
            continue
        timeDirs.append(path)
timeDirs = [ x[1] for x in sorted(zip(times,timeDirs)) ]
times.sort()
Ntimes = len(times)

# read the mesh (cell centered points)
xcc = np.zeros(N) # these are the cell-centered coordinates
ycc = np.zeros(N)
zcc = np.zeros(N)
meshfile = os.path.join(timeDirs[0],prefix+'.mesh')
if verbose: print 'Reading mesh from',meshfile
with open(meshfile,'r') as f:
    for _ in range(8): f.readline()
    assert(int(f.readline()) == N)
    for i in range(N):
        xcc[i] = float(f.readline())
    for i in range(N):
        ycc[i] = float(f.readline())
    for i in range(N):
        zcc[i] = float(f.readline())
xcc = xcc.reshape((Nx,Ny,Nz),order='F')
ycc = ycc.reshape((Nx,Ny,Nz),order='F')
zcc = zcc.reshape((Nx,Ny,Nz),order='F')
if verbose:
    print '  range x: ',[np.min(xcc),np.max(xcc)]
    print '  range y: ',[np.min(ycc),np.max(ycc)]
    print '  range z: ',[np.min(zcc),np.max(zcc)]

# create the actual mesh
if Nx > 1:
    dx = xcc[1,0,0] - xcc[0,0,0]
    if verbose: print 'dx=',dx
    #x1 = xcc[0,0,0]-dx/2 + np.arange(Npx)*dx
    x1 = xcc[:,0,0]
else:
    dx = 1.0
    x1 = [xcc[0,0,0]]
if Ny > 1:
    dy = ycc[0,1,0] - ycc[0,0,0]
    if verbose: print 'dy=',dy
    #y1 = ycc[0,0,0]-dy/2 + np.arange(Npy)*dy
    y1 = ycc[0,:,0]
else:
    dy = 1.0
    y1 = [ycc[0,0,0]]
if Nz > 1:
    dz = zcc[0,0,1] - zcc[0,0,0]
    if verbose: print 'dz=',dz
    #z1 = zcc[0,0,0]-dz/2 + np.arange(Npz)*dz
    z1 = zcc[0,0,:]
else:
    dz = 1.0
    z1 = [zcc[0,0,0]]
X,Y,Z = np.meshgrid(x1,y1,z1,indexing='ij')

# process all velocity fields
outdir = os.path.join(samplingDir,'VTK')
if not os.path.isdir(outdir):
    os.makedirs(outdir)
laststr = ''
for itime,tdir in enumerate(timeDirs):
    solnfile = os.path.join(tdir,prefix+'.000.U')
    newstr = '[{:2d}%] Processing {:s}'.format(int(100.0*itime/Ntimes),solnfile)
    sys.stdout.write('\r'+len(laststr)*' ')
    sys.stdout.write(newstr)
    # read solution file
    U = np.zeros(N)
    V = np.zeros(N)
    W = np.zeros(N)
    with open(solnfile,'r') as f:
        for _ in range(4): f.readline()
        lines = f.readlines()
    data = [ float(val) for val in lines ]
    U[:] = data[:N]
    V[:] = data[N:2*N]
    W[:] = data[2*N:3*N]
    U = U.reshape((Nx,Ny,Nz),order='F')
    V = V.reshape((Nx,Ny,Nz),order='F')
    W = W.reshape((Nx,Ny,Nz),order='F')
    # write binary VTK
    outfile = os.path.join(outdir,prefix+'_'+str(itime)+'.vtk')
    VTKwriter.vtk_write_structured_points( open(outfile,'wb'), #binary mode
        Nx, Ny, Nz,
        [ U,V,W ],
        datatype=['vector'],
        dx=dx, dy=dy, dz=dz,
        dataname=['U'],
        origin=[x1[0],y1[0],z1[0]],
        indexorder='ijk')

sys.stdout.write('\n...done\n')
