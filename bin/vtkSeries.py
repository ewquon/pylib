#!/usr/bin/env python
#
# for creating symlinks to output from the OpenFOAM sample utility
# expected directory structure:
#   postProcessing/surfaces/<time>/<var>_<setName>.vtk
# result:
#   postProcessing/surfaces/<setName>/<var>_<time>.vtk
#
import os

dirlist = []
timesteps = []

#dirs=*.*
#curdir=`pwd`
dirs = os.walk('.').next()[1]
for d in dirs:
    try: 
        step = float(d) # need this to verify this is a time-step dir!
        dirlist.append(d)
        timesteps.append(step)
    except ValueError: pass

extMapping = dict()
extMapping['xy'] = 'xyz'

def tname(tval):
    # note: paraview doesn't seem to handle floats well...
    #return '%d' % (tval*10)
    #return '%d' % (tval*1000)
    return '%d' % (tval)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

sampleNames = []
varNames = []
extNames = []
for timestep_dir in dirlist:
    print 'Processing', timestep_dir
    for f in [ f for f in os.listdir(timestep_dir) if os.path.isfile(os.path.join(timestep_dir,f)) ]:
        if f.startswith('.'):
            continue
#        fsplit = f.split('.')
#        fbasename = '.'.join(fsplit[0:-1])
#        ext = fsplit[-1]
        fbasename,ext = os.path.splitext(f)
        ext = ext[1:]
        fbasesplit = fbasename.split('_')
        var = fbasesplit[0]
        name = '_'.join(fbasesplit[1:])
        print ' ',f,'( name=',name,' var=',var,' ext=',ext,')'
        if name=='':
            name = 'timeSeries'
        if not name in sampleNames:
            sampleNames.append(name)
            if not os.path.exists(name): os.makedirs(name)
        if not var in varNames:
            varNames.append(var)
        if not ext in extNames:
            extNames.append(ext)

if not len(extNames)==1:
    print 'Don''t know how to handle different extensions',extNames
if ext in extMapping:
    extNew = extMapping[ext]
else: extNew = ext

indices = sorted(range(len(timesteps)), key=lambda k: timesteps[k])
for sample in sampleNames:
    for var in varNames:
        for i in range(len(timesteps)):
            idx = indices[i]
            dname = dirlist[idx]#.split()
            if sample=='timeSeries':
                src = os.path.join( os.getcwd(), dname, var+'.'+ext )
                #dest = sample + os.sep + '%s_%s.%s' % (var,tname(timesteps[idx]),extNew)
                dest = sample + os.sep + '%s_%s.%s' % (var,i,extNew)
            else:
                #src = os.path.join( os.getcwd(), dname, sample+'_'+var+'.'+ext )
                # result, e.g.: /Users/equon/wndpltdesign/inflow/5mps_ab/postProcessing/surfaces/9000/slice_cross_stream_U.vtk
                src = os.path.join( os.getcwd(), dname, var+'_'+sample+'.'+ext )
                #dest = sample + os.sep + '%s_%s.%s' % (var,tname(timesteps[idx]),extNew)
                dest = sample + os.sep + '%s_%s.%s' % (var,i,extNew)
            print dest,'-->',src
            try:
                os.symlink(src,dest)
            except OSError: pass

