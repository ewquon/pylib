import numpy as np
import sys

class my_vtk:

    verbose = True

    def __init__(self,fname=''):

        self.filename = fname

        if not fname=='':
            self.read_vtk_surf(fname)


    ##
    ## Read polygonal surface data
    ##
    def read_vtk_surf(self,fname):

        with open(fname,'r') as f:

            f.readline()
            self.title = f.readline().strip()
            line = f.readline()
            self.datatype = line.strip().lower()
            if not self.datatype=='ascii':
                print 'Data format',self.datatype,'not handled yet!'
                return None

            line = f.readline().split()
            assert( line[0].lower()=='dataset' )
            self.geomtype = line[1].lower()

            if self.geomtype=='polydata':

                line = f.readline().split()
                assert( line[0].lower()=='points' )
                assert( line[2].lower()=='float' )
                self.npts = int(line[1])

                self.pts = np.zeros((self.npts,3))
                for i in range(self.npts):
                    self.pts[i,:] = [ float(val) for val in f.readline().split() ]

                assert( f.readline().strip()=='' )

                line = f.readline().split()
                assert( line[0].lower()=='polygons' )
                self.npoly = int(line[1])
                nsize = int(line[2])
                self.npolypts = nsize/self.npoly - 1
                assert( nsize%self.npoly==0 )

                self.connectivity = np.zeros((self.npoly,self.npolypts),dtype=int)
                self.polycenters = None
                for i in range(self.npoly):
                    self.connectivity[i,:] = [ int(val) for val in f.readline().split()[1:] ]

                line = f.readline().split()
                assert( line[0].lower()=='point_data' )
                line = f.readline().split()
                assert( line[0].lower()=='field' )
                line = f.readline().split()
                assert( line[0].lower()=='p' )
                assert( int(line[2])==self.npts )
                assert( line[3].lower()=='float' )

                pdata = []
                for line in f:
                    pdata += line.split()
                self.data = np.array([ float(p) for p in pdata ])
                assert( len(self.data)==self.npts )

            else:
                print 'Geometry type',self.geomtype,'not handled yet!'
                return None
            
    ##
    ## Return cell-centered coordinates, calculating if necessary
    ##
    def cc(self):

        if not self.polycenters:
            #if self.verbose: print "VTK: calculating cell centers"
            if self.verbose: sys.stdout.write(' [calculating cell centers]')
            self.x_cc = np.zeros((self.npoly,3))
            self.data_cc = np.zeros(self.npoly)

            for i in range(self.npoly):
                conn = self.connectivity[i,:]
                for j in range(3):
                    self.x_cc[i,j] = np.mean(self.pts[conn,j])
                self.data_cc[i] = np.mean(self.data[conn])

            self.polycenters = True # don't need to recalculate this?

        return self.x_cc[:,0], self.x_cc[:,1], self.x_cc[:,2], self.data_cc

    # return cell-centered coordinates of cells that straddle the specified location
    def cc_at_const(self,slice_loc,slice_dir):
        x,y,z,f = self.cc()
        idx = self.faces_at_const(slice_loc,slice_dir)
        xerr = self.x_cc[idx,0] - slice_loc
        yerr = self.x_cc[idx,1] - slice_loc
        zerr = self.x_cc[idx,2] - slice_loc
        errors = [ np.std(xerr), np.std(yerr), np.std(zerr) ]
       # if self.verbose: sys.stdout.write(' [pos err=%g]'%errors[slice_dir])
        if self.verbose: sys.stdout.write(' [pos err=%g]\n' % errors[slice_dir])
        return x[idx],y[idx],z[idx],f[idx]

    def cc_at_x(self,x): return self.cc_at_const(x,0)
    def cc_at_y(self,y): return self.cc_at_const(y,1)
    def cc_at_z(self,z): return self.cc_at_const(z,2)

    # return IDs of faces (assuming 2D polygons)
    def faces_at_const(self,slice_loc,slice_dir):
        idx = []
        for i in range(self.npoly):
            conn = self.connectivity[i,:]
            minpos = np.min(self.pts[conn,slice_dir])
            maxpos = np.max(self.pts[conn,slice_dir])
            if slice_loc >= minpos and slice_loc < maxpos: idx.append(i)
        return idx

    def faces_at_x(self,x): return self.faces_at_const(x,0)
    def faces_at_y(self,y): return self.faces_at_const(y,1)
    def faces_at_z(self,z): return self.faces_at_const(z,2)
    

    ##
    ## Calculate face normals
    ##
    def face_normals(self,cellIDs=''):
        if len(cellIDs)==0: 
            cells = self.connectivity[:,:]
            Ncells = self.npoly
        else:
            cells = self.connectivity[cellIDs,:]
            Ncells = len(cellIDs)

        Sf = np.zeros((Ncells,3))
        magSf = np.zeros(Ncells)

        def norm(v): return np.dot(v,v)**0.5

        for icell in range(Ncells):
            ptIDs = cells[icell,:]
            vecA = self.pts[ptIDs[1],:] - self.pts[ptIDs[0],:]
            vecB = self.pts[ptIDs[2],:] - self.pts[ptIDs[1],:]
            Sf[icell,:] = np.cross(vecA,vecB)
            magSf[icell] = norm(Sf[icell,:]) # twice the area spanned by the vectors
            Sf[icell,:] /= magSf[icell] # unit normal

            if len(ptIDs) > 3: # quad cell, assuming no higher-order polygons
                vecC = self.pts[ptIDs[3],:] - self.pts[ptIDs[2],:]
                vecD = self.pts[ptIDs[0],:] - self.pts[ptIDs[3],:]
                Sf2 = np.cross(vecC,vecD)
                mag2 = norm(Sf2) # twice the area spanned by the vectors
                Sf2 /= mag2 # unit normal
                #assert(all(np.abs(Sf2-Sf[i,:]) < 1e-6)) # make sure both vectors point in the same direction
                Sf[icell,:] = 0.5*( Sf[icell,:] + Sf2 ) # use average of direction
                magSf[icell] += mag2 # add area of second half of polygon

            magSf[icell] /= 2 # get correct area

            # assuming all cells have same winding and are pointing outward...
            # TODO: check this!

        return Sf, magSf


