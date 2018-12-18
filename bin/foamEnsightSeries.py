#!/usr/bin/env python
#
# Create symbolic links for specified directories to form a single Ensight .case
# file with the complete time series
#
from __future__ import print_function
import sys
import os
import glob
import numpy as np

# HARD CODE INPUTS FOR NOW
prefixMeshFile = '_U.mesh'
prefixSolnFile = '_U.000.U'
prefixSolnStr  = '_U.*****.U'
prefixNewFile  = '_U.{:05d}.U'

if len(sys.argv) > 1:
    postdir = sys.argv[1].rstrip(os.sep)
else:
    sys.exit('Specify array postprocessing directory with structure: arraySampleDir/timeDir/*.case')
if not os.path.isdir(postdir):
    sys.exit('Not a valid directory: '+postdir)

if len(sys.argv) > 2:
    prefix = sys.argv[2]
else:
    prefix = os.path.split(postdir)[1]
meshFile = prefix + prefixMeshFile
solnFile = prefix + prefixSolnFile
solnStr  = prefix + prefixSolnStr
newFile  = prefix + prefixNewFile

outdir = postdir + '_series'
if os.path.isdir(outdir):
    print('Warning: time series directory',outdir,'already exists')
else:
    os.makedirs(outdir)

# read available time directories (assume equally spaced)
timeList = []
tdirList = []
for d in glob.glob(os.path.join(postdir,'*')):
    if os.path.isdir(d):
        tdir = os.path.split(d)[-1]
        try:
            tval = float(tdir)
        except ValueError:
            continue
        timeList.append(tval)
        tdirList.append(os.path.abspath(d))

timeArray = np.array(timeList)
timeSortOrder = np.argsort(timeArray)
timeArray = timeArray[timeSortOrder]
tdirList = [ tdirList[i] for i in timeSortOrder ]
#for t,tdir in zip(timeArray,tdirList):
#    print('t=',t,':',tdir)
Ntimes = len(tdirList)

timeDiff = np.diff(timeArray)
if not np.min(timeDiff) == np.max(timeDiff):
    print('Warning: sampling intervals are variable?')

# create symlinks
havemesh = False
for itime in range(Ntimes):
    if not havemesh:
        os.symlink(
                os.path.join(tdirList[itime],meshFile),
                os.path.join(outdir,meshFile)
                )
        havemesh = True
    srcfile = os.path.join(tdirList[itime],solnFile)
    if not os.path.isfile(srcfile):
        raise IOError('Source file not found: '+srcfile)
    tgtfile = os.path.join(outdir,newFile.format(itime))
    print(tgtfile,'-->',srcfile)
    os.symlink(srcfile,tgtfile)

# create case file for time series
caseTemplate="""FORMAT
type: ensight gold

GEOMETRY
model:        1     {meshFile:s}

VARIABLE
vector per node:            1       U               {solnStr:s}

TIME
time set:                      1
number of steps:               {Ntimes:d}
filename start number:         0
filename increment:            1
time values:
"""
targetMesh = os.path.join(outdir,meshFile)
targetSoln = os.path.join(outdir,solnStr)
with open(outdir+'.case','w') as f:
    f.write(
            caseTemplate.format(meshFile=targetMesh,solnStr=targetSoln,Ntimes=Ntimes)
            )
    for t in timeArray:
        f.write('{:.5g}\n'.format(t))


