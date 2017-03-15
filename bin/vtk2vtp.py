#!/usr/bin/env python
# based on https://gist.github.com/thomasballinger/1281457
import sys
import os

import vtk

def vtk2vtp(inputFile, outputFile, verbose=True):
    """Convert legacy VTK file to new VTK Polygonal Data format

    Note: binary files written on one system may not be readable on
    another.
    """
    if not os.path.isfile(inputFile):
        raise IOError('Input file not found: '+inputFile)

    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(inputFile)
    reader.Update()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outputFile)
    #writer.SetFileTypeToBinary()
    writer.SetInputData(reader.GetOutput())
    writer.Update()

    if verbose:
        print 'Wrote',outputFile

#=======================================================================
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print 'specify file(s) to batch convert'
        sys.exit()

    for vtkfile in sys.argv[1:]:
        if not vtkfile[-4:] == '.vtk':
            print 'Skipping',vtkfile
            continue

        vtk2vtp(vtkfile, vtkfile[:-4]+'.vtp')

