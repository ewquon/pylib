#!/usr/bin/env python
import sys
import os

import vtk

def combineVtkCellData(*args, **kwargs):
    """Combine multiple legacy VTK files into a single file with the VTK
    Polygonal Data format. Files are assumed to contain only cell data
    arrays.
    """
    outputFile = kwargs.get('outputFile','new.vtp')
    verbose = kwargs.get('verbose',True)
    data = None

    for inputFile in args:
        if not os.path.isfile(inputFile):
            raise IOError('Input file not found: '+inputFile)
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(inputFile)
        reader.Update()

        if data is None:
            data = vtk.vtkPolyData()
            data.DeepCopy(reader.GetOutput())
            celldata = data.GetCellData()
        else:
            appendcelldata = reader.GetOutput().GetCellData()
            for Iarray in range(appendcelldata.GetNumberOfArrays()):
                newarray = appendcelldata.GetAbstractArray(Iarray)
                celldata.AddArray(newarray)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outputFile)
    #writer.SetFileTypeToBinary()
    writer.SetInputData(data)
    writer.Update()

    if verbose:
        print 'Wrote',outputFile

#=======================================================================
if __name__ == '__main__':
    if len(sys.argv) <= 2:
        print 'specify two or more files to combine'
        sys.exit()

    combineVtkCellData(*sys.argv[1:], outputFile='combined.vtp', verbose=False)

