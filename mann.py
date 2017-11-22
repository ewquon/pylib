#!/usr/bin/env python
#
# Mann turbulent field data processing module (output from WaSP/windsimu)
# written by Eliot Quon (eliot.quon@nrel.gov) - 2017-07-06
#
import sys,os
import numpy as np

import binario
import VTKwriter

class bin:

    def __init__(self,
            fname,
            Nx=1,Ny=1,Nz=1,Ncomp=1,
            Lx=None,Ly=None,Lz=None,
            verbose=True):

        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.Ncomp = Ncomp

        if Lx is None: Lx = float(Nx)
        if Ly is None: Ly = float(Ny)
        if Lz is None: Lz = float(Nz)
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dx = Lx/Nx
        self.dy = Ly/Ny
        self.dz = Lz/Nz

        if verbose:
            print 'Domain extents: ',[self.Lx,self.Ly,self.Lz]
            print 'Cell spacings: ', [self.dx,self.dy,self.dz]

        self._readBinary(fname)

        if verbose:
            print 'Velocity component ranges:'
            for i in range(self.Ncomp):
                print '  u'+str(i)+': ',[np.min(self.u[:,:,:,i]),np.max(self.u[:,:,:,i])]

    def _readBinary(self,fname):
        N = self.Nx * self.Ny * self.Nz
        self.u = np.zeros((self.Nx,self.Ny,self.Nz,self.Ncomp))
        with binario.binaryfile(fname) as f:
            def readField():
                data = f.read_real4(N)
                return np.array(data).reshape((self.Nx,self.Ny,self.Nz),order='C')
            for i in range(self.Ncomp):
                self.u[:,:,:,i] = f.readField()

