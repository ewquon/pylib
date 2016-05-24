#!/usr/local/bin/python
import numpy as np

class domain:
    def __init__(self,Lx,Ly,Lz,dsmax=1.0,dt=1.0):
        self.Ndomain = 0
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dt = dt

        self.name = []
        self.level = []
        self.cellvol = []
        self.Ncells0 = []
        self.Nx = []
        self.Ny = []
        self.Nz = []
        self.dx = []
        self.dy = []
        self.dz = []
        self.ds = []
        self.volume0 = [] # volume of each input domain
        self.volume = [] # volume in each level
        
        self.add_subdomain(level=0,Xref=Lx,Yref=Ly,Zref=Lz,ds=dsmax,name='full domain')

    def __repr__(self):
        dsmax = np.log10( np.max((self.dx,self.dy,self.dz)) )
        dsmin = np.log10( np.min((self.dx,self.dy,self.dz)) )
        output = 'Input configuration:'
        for i,domname in enumerate(self.name):
            output += '\n  Level {} : dx,dy,dz={} {} {} ncells={} vol={}  {}'.format( \
                    self.level[i], \
                    self.dx[i], self.dy[i], self.dz[i], \
                    self.Ncells0[i], \
                    self.volume0[i], \
                    domname )
        output += '\nmax/min/diff in scale: {} {} {}'.format(dsmax,dsmin,dsmax-dsmin)
        return output

    def add_subdomain(self,level=0,Xref=1.0,Yref=1.0,Zref=1.0,\
            ds=-1,dx=-1,dy=-1,dz=-1,\
            name='subdomain'):
        self.name.append( name )
        self.level.append( level )

        if ds > 0:
            dx = ds
            dy = ds
            dz = ds
        else:
            assert( dx > 0 and dy > 0 and dz > 0 )
        #self.ds.append( ds )
        self.dx.append( dx )
        self.dy.append( dy )
        self.dz.append( dz )
        #self.cellvol.append( ds**3 )
        self.cellvol.append( dx*dy*dz )
        self.ds.append( self.cellvol[-1]**(1./3.) )
        self.Nx.append( int(Xref / dx) )
        self.Ny.append( int(Yref / dy) )
        self.Nz.append( int(Zref / dz) )

        #print 'adding subdomain',self.Ndomain,'L|W|H =',Xref,Yref,Zref,\
        #        'Nx|y|z =',self.Nx[-1],self.Ny[-1],self.Nz[-1]
        self.Ncells0.append( self.Nx[-1]*self.Ny[-1]*self.Nz[-1] )
        self.volume0.append( self.Ncells0[-1] * self.cellvol[-1] )
        self.Ndomain += 1

    def check_total_volume(self):
        print 'checking total volume...'
        maxlevel = np.max(self.level)
        print '  max level =',maxlevel
        self.volume = np.zeros(maxlevel+1)
        self.Ncells = np.array(self.Ncells0)

        # calculate total volume per level
        for i,vol in enumerate(self.volume0):
            self.volume[self.level[i]] = vol

        # adjust volume and Ncells so we don't double count
        for i in range(1,self.Ndomain):
            ilevel = self.level[i]
            self.volume[ilevel-1] -= self.volume[ilevel]
            self.Ncells[ilevel-1] -= int(self.volume[ilevel]/self.cellvol[ilevel-1])

        vol_calc = np.sum(self.volume)
        print '  volumes per domain :',self.volume
        print '  total volume = {:g} {:g} (err={:g}%)'.format( \
                self.volume0[0], vol_calc, 100*(self.volume0[0]-vol_calc)/self.volume0[0] )
        print '  cells per domain :',self.Ncells
        print '  total cells = {:g}'.format( np.sum(self.Ncells) )

D = 100.0   # ref length = rotor diameter
cref = 1.0  # ref chord length
tref = 0.1  # ref airfoil thickness
dom = domain(5000,5000,1000,dsmax=10.0,dt=1e-3)
dom.add_subdomain( level=1, Xref=5.*D , Yref=D   , Zref=D  , ds=1.25,  name='far wake')
dom.add_subdomain( level=2, Xref= D/2., Yref=D   , Zref=D  , ds=0.05,  name='near wake')
#dom.add_subdomain( level=3, Xref= cref, Yref=tref, Zref=D/2, ds=5.e-6, name='blade 1')
#dom.add_subdomain( level=3, Xref= cref, Yref=tref, Zref=D/2, ds=5.e-6, name='blade 2')
#dom.add_subdomain( level=3, Xref= cref, Yref=tref, Zref=D/2, ds=5.e-6, name='blade 3')
dom.add_subdomain( level=3, Xref= cref, Yref=tref, Zref=D/2, dx=0.01, dy=0.01, dz=5.e-6, name='blade 1')
dom.add_subdomain( level=3, Xref= cref, Yref=tref, Zref=D/2, dx=0.01, dy=0.01, dz=5.e-6, name='blade 2')
dom.add_subdomain( level=3, Xref= cref, Yref=tref, Zref=D/2, dx=0.01, dy=0.01, dz=5.e-6, name='blade 3')
print dom

dom.check_total_volume()
