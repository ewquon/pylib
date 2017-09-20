#!/usr/bin/env python
#
# Script for quickly checking rotor loads output from SOWFA
#
# Eliot Quon <eliot.quon@nrel.gov>
#
import numpy as np
import matplotlib.pyplot as plt
from SOWFA.postProcessing import turbineOutput

Nturb = 4
toffset = 21600 #None
Omega = 10.0 # RPM
dt = 0.04

#===============================================================================

Trotate = 1.0/(Omega/60.0) # s / cycle
plotskip = int(Trotate/4.0/dt) # 4 overlaid snapshots per nominal rotor period

rotorOutputs = turbineOutput.readRotorOutputs(Nturb=Nturb,toffset=toffset)
#bladeOutputs = turbineOutput.readBladeOutputs(Nturb=Nturb,toffset=toffset)
#towerOutputs = turbineOutput.readTowerOutputs(Nturb=Nturb,toffset=toffset)

def plotProperty(dflist,propertyName,turblist=range(Nturb)):
    """Plot bulk (rotor-level) quantities versus time"""
    selectedProperty = None
    for outputName in list(dflist[0]):
        if propertyName.lower() in outputName.lower():
            selectedProperty = outputName
            break
    if selectedProperty is not None:
        print 'Plotting',selectedProperty
        plt.figure()
        for iturb in turblist:
            dflist[iturb][selectedProperty].plot(label='turbine {:d}'.format(iturb))
        plt.ylabel(selectedProperty)
        plt.title(propertyName)
        plt.legend(loc='best',fontsize='small')
    else:
        print 'Property "{:s}" not found'.format(propertyName)

def plotDistributedProperty(dflist,propertyName,turblist=range(Nturb),plotskip=plotskip):
    """Plot distributed quantities over the turbine component span

    The unsteady variation in a quantity is lightly plotted, while the time
    average is plotted in bold color.
    """
    selectedColumns = []
    for outputName in list(dflist[0]):
        if propertyName.lower() in outputName.lower():
            selectedColumns.append(outputName)
    numColumns = len(selectedColumns)
    if numColumns > 0:
        outputName = '_'.join(selectedColumns[0].split('_')[:-1])
        print 'Plotting {:s} ({:d} points)'.format(outputName,numColumns)
        x = np.arange(numColumns)

        plt.figure()
        for iturb in turblist:
            values = dflist[iturb].iloc[0].ix[selectedColumns]
            line, = plt.plot(x,values,label='turbine {:d}'.format(iturb),alpha=0.1)
            color = line._color
            for irow in range(1,len(dflist[iturb].index),plotskip):
                values = dflist[iturb].iloc[irow].ix[selectedColumns].values
                plt.plot(x,values,color=color,alpha=0.1)
        plt.ylabel(outputName)
        plt.title(propertyName)
        plt.legend(loc='best',fontsize='small')
    else:
        print 'Distributed property "{:s}" not found'.format(propertyName)

print 'rotor outputs :',list(rotorOutputs[0])
#print 'blade outputs :',list(bladeOutputs[0])
#print 'tower outputs :',list(towerOutputs[0])

plotProperty(rotorOutputs,'Power')
plotProperty(rotorOutputs,'Rotation Rate')
plotProperty(rotorOutputs,'Rotor Axial Force')

#plotDistributedProperty(bladeOutputs,'Axial Force')

plt.show()

