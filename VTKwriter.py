#!/usr/local/bin/python
#
# VTK I/O helper module
# written by Eliot Quon (eliot.quon@nrel.gov)
#
import sys
import numpy as np

# constant
vtk_header = '# vtk DataFile Version 2.0'
vtk_datasettype = 'structured_points'
vtk_datatype = 'float'

# input
vtk_description = 'really cool data'
#vtk_format = 'ascii'

"""NOT USED:"""
"""vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"""
def vtk_write_rectilinear_grid(f,x,y,z,u,v,w,dataname='velocity'):# {{{
    """ write a dataset with regular topology to file f
    x,y,z may be 1-D arrays or full meshgrids
    Inputs are written with x increasing fastest, then y, then z
    3D vector data are assumed
    NOT TESTED YET
    """
    f.write(vtk_header+'\n')
    f.write(vtk_description+'\n')
    f.write(vtk_format.upper()+'\n')
    if len(x.shape) > 1:
        nx,ny,nz = x.shape
        x = x[:,0,0]
        y = y[0,:,0]
        z = z[0,0,:]
    else:
        nx = len(x)
        ny = len(y)
        nz = len(z)
    f.write('DATASET RECTILINEAR_GRID')
    f.write('DIMENSIONS {:d} {:d} {:d}\n'.format(nx,ny,nz))

    f.write('X_COORDINATES {:d} {:s}\n'.format(nx,vtk_datatype))
    for xi in x:
        f.write(' {:f}'.format(xi))
    f.write('\n')
    f.write('Y_COORDINATES {:d} {:s}\n'.format(ny,vtk_datatype))
    for yi in y:
        f.write(' {:f}'.format(yi))
    f.write('\n')
    f.write('Z_COORDINATES {:d} {:s}\n'.format(nz,vtk_datatype))
    for zi in z:
        f.write(' {:f}'.format(zi))
    f.write('\n')

    f.write('POINT_DATA {:d}\n'.format(nx*ny*nz))
    f.write('VECTORS {:s} {:s}\n'.format(dataname,vtk_datatype))
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                f.write(' {:f} {:f} {:f}\n'.format(u[j,i,k],v[j,i,k],w[j,i,k]))# }}}
"""^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"""

def vtk_write_structured_points(f,nx,ny,nz, # {{{
        data, datatype=['vector'],
        ds=None,dx=None,dy=None,dz=None,
        origin=(0.0,0.0,0.0),
        dataname=[],
        indexorder='jik'):
    """ write a dataset with regular topology to file f
    'ds' is the default grid spacing; 'dx','dy','dz' may be specified to override.
    Inputs are written with x increasing fastest, then y, then z.
      Acceptable types: 'vector','scalar'
    Additional vectors may be specified as additional scalar arguments; example:
      vtk_write_structured_points(f,nx,ny,nz,[u,v,w,u1,v2,w2],ds=1.0,dataname=['mean','fluctuation'])
    """
    # caluclate grid spacings if needed
    if ds:
        if not dx: dx = ds
        if not dy: dy = ds
        if not dz: dz = ds
    else:
        assert( dx > 0 and dy > 0 and dz > 0 ) 

    # replace shorthand names
    if type(dataname)==str: dataname = [dataname]
    Nvector = 0
    Nscalar = 0
    Nvalues = 0
    for i,name in enumerate(datatype):
        if name[0].lower() == 'v':
            datatype[i] = 'vector'
            Nvector += 1
            Nvalues += 3
        elif name[0].lower() == 's':
            datatype[i] = 'scalar'
            Nscalar += 1
            Nvalues += 1
        else:
            print 'unrecognized data type',name

    # sanity checks
    #Nscalars = len(data)
    #Nvecs = Nscalars/3
    #assert( Nvecs*3 == Nscalars )
    assert( len(data) == Nvalues )

    # write header
    f.write(vtk_header+'\n')
    f.write(vtk_description+'\n')
    if 'b' in f.mode:
        binary = True
        import struct
        f.write('BINARY\n')
    else:
        binary = False
        f.write('ASCII\n')
    f.write('DATASET STRUCTURED_POINTS\n')

    # write out mesh descriptors
    f.write('DIMENSIONS {:d} {:d} {:d}\n'.format(nx,ny,nz))
    f.write('ORIGIN {:f} {:f} {:f}\n'.format(origin[0],origin[1],origin[2]))
    f.write('SPACING {:f} {:f} {:f}\n'.format(dx,dy,dz))

    # write out data
    f.write('POINT_DATA {:d}\n'.format(nx*ny*nz))
    idx = 0 # data list index
    for idata,outputtype in enumerate(datatype):

        if outputtype=='vector':
            u,v,w = data[idx], data[idx+1], data[idx+2]
            idx += 3
        elif outputtype=='scalar':
            u = data[idx]
            idx += 1
        else: continue

        try: name = dataname[idata]
        except IndexError: name = outputtype+str(idata)

        if outputtype=='vector':
            f.write('{:s}S {:s} {:s}\n'.format(outputtype.upper(),name,vtk_datatype))
            if indexorder=='jik': # TTUDD indexing convention:
                if binary:
                    for k in range(nz):
                        for j in range(ny):
                            for i in range(nx):
                                #f.write(struct.pack('fff',u[j,i,k],v[j,i,k],w[j,i,k])) # native endianness
                                f.write(struct.pack('>fff',u[j,i,k],v[j,i,k],w[j,i,k])) # big endian
                else: #ascii
                    for k in range(nz):
                        for j in range(ny):
                            for i in range(nx):
                                #f.write(' {:f} {:f} {:f}\n'.format(u[i,j,k],v[i,j,k],w[i,j,k]))
                                f.write(' {:f} {:f} {:f}\n'.format(u[j,i,k],v[j,i,k],w[j,i,k]))
            else: # assume index order is i,j,k
                if binary:
                    for k in range(nz):
                        for j in range(ny):
                            for i in range(nx):
                                #f.write(struct.pack('fff',u[i,j,k],v[i,j,k],w[i,j,k])) # native endianness
                                f.write(struct.pack('>fff',u[i,j,k],v[i,j,k],w[i,j,k])) # big endian
                else: #ascii
                    for k in range(nz):
                        for j in range(ny):
                            for i in range(nx):
                                #f.write(' {:f} {:f} {:f}\n'.format(u[i,j,k],v[i,j,k],w[i,j,k]))
                                f.write(' {:f} {:f} {:f}\n'.format(u[i,j,k],v[i,j,k],w[i,j,k]))

        elif outputtype=='scalar':
            f.write('{:s}S {:s} {:s}\n'.format(outputtype.upper(),name,vtk_datatype))
            f.write('LOOKUP_TABLE default\n')
            if binary:
                for k in range(nz):
                    for j in range(ny):
                        for i in range(nx):
                            #f.write(struct.pack('f',u[j,i,k])) # native endianness
                            f.write(struct.pack('>f',u[j,i,k])) # big endian
            else:
                for k in range(nz):
                    for j in range(ny):
                        for i in range(nx):
                            f.write(' {:f}\n'.format(u[j,i,k]))


# }}}

