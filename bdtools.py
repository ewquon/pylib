#!/usr/local/bin/python
import plotcols
import matplotlib.pyplot as plt
import numpy as np

colMap = dict()
colMap['Fx'] = 1
colMap['Fy'] = 2
colMap['Fz'] = 3
colMap['Mx'] = 4
colMap['My'] = 5
colMap['Mz'] = 6
colMap['x'] = 1
colMap['y'] = 2
colMap['z'] = 3
colMap['p'] = 4
colMap['q'] = 5
colMap['r'] = 6

GLLpts = []
beamnodes = []

#---------------------------------------------------------------
# wrappers around plotcols

# in the IEC ref frame
def plotLoad(varnames,fname='load.out'):
    cols = []
    for name in varnames:
        loadType = name[0]
        loadDir  = name[1]
        beamSegment = int(name[2:])
        cols.append( 6*(beamSegment-1) + colMap[loadType+loadDir] )
    plotcols.plot(fname,cols,varnames,xname='t',title=fname)

# in the IEC ref frame
def plotDisp(varnames,fname='rel_disp.out'):
    cols = []
    for name in varnames:
        dispType = name[0]
        node = int(name[1:])
        cols.append( 6*(node-1) + colMap[dispType] )
    plotcols.plot(fname,cols,varnames,xname='t',title=fname)

#---------------------------------------------------------------
# helper functions

def read_beam_input(searchStr,thenSkip=0,fname='beam.input.sum'):
    with open(fname,'r') as f:
        line = f.readline()
        while not line.startswith(searchStr):
            line = f.readline()
        for i in range(thenSkip): f.readline()
        line = f.readline().split()
        nCols = len(line)
        data = []
        for col in range(nCols): data.append( [ float(line[col]) ] )
        line = f.readline().split()
        while not len(line)==0:
            for col in range(nCols): data[col].append( float(line[col]) )
            line = f.readline().split()
    return data


#---------------------------------------------------------------
# section plots

def plotSectionLoads(varnames,sections=[],fname='load.out'):
    global GLLpts
    if len(GLLpts)==0: GLLpts = read_beam_input('Quadrature point')[1]

    if len(sections)==0: sections = range(len(GLLpts))

    for varname in varnames:
        cols = []
        for sect in sections:
            cols.append( 6*(sect-1) + colMap[varname] )

        #plt.figure()
        with open(fname,'r') as f:
            for line in f: pass # plot latest only
            line = line.split()
            plt.plot([ GLLpts[i] for i in sections ], \
                     [ line[i] for i in cols ], \
                     label=varname
                     #color='0.5'
                    )

    plt.xlabel('x')
    #plt.ylabel(varname)
    plt.title(fname)
    plt.legend(loc='best')
    plt.show(block=False)


def plotSectionDisp(varnames,nodes=[],fname='rel_disp.out'):
    global beamnodes
    if len(beamnodes)==0: beamnodes = read_beam_input('Initial pos',thenSkip=1)[1]

    if len(nodes)==0: nodes = range(len(beamnodes))

    for varname in varnames:
        cols = []
        for node in nodes:
            cols.append( 6*(node-1) + colMap[varname] )

        #plt.figure()
        with open(fname,'r') as f:
            for line in f: pass # plot latest only
            line = line.split()
            plt.plot([ beamnodes[i] for i in nodes ], \
                     [ line[i] for i in cols ], \
                     label=varname
                     #color='0.5'
                    )

    plt.xlabel('x')
    #plt.ylabel(varname)
    plt.title(fname)
    plt.legend(loc='best')
    plt.show(block=False)

