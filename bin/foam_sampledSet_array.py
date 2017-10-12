#!/usr/bin/env python
import numpy as np

def sampleCoords(origin,pointsDensity,spanBox):
    nTotalSamples = np.prod(pointsDensity[:3])

    sampleCoords = np.zeros((nTotalSamples,3))

    deltax = spanBox[0] / (pointsDensity[0] + 1)
    deltay = spanBox[1] / (pointsDensity[1] + 1)
    deltaz = spanBox[2] / (pointsDensity[2] + 1)

    p = 0
    for k in range(1,pointsDensity[2]+1):
        for j in range(1,pointsDensity[1]+1):
            for i in range(1,pointsDensity[0]+1):
                t = np.array([deltax*i, deltay*j, deltaz*k])
                sampleCoords[p,:] = origin + t
                p += 1

    for i in range(nTotalSamples):
        print sampleCoords[i,:]

#-------------------------------------------------------------------------------

origin = [495.0,1415.0,5.0]
pointsDensity = [1,16,16]
spanBox = [20,170,170]

sampleCoords(origin,pointsDensity,spanBox)
 
