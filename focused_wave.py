#!/usr/bin/python
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

class focused_wave:

    #
    # defaults
    #
    x_wec = 0.0 # device location
    x_inlet = -600.0 # (near) inlet
    dt = 0.02
    surfDir = 'waterSurface'
    MLERinput = 'coeffs.txt'
    toffset = 150.0 # time after which the focused wave reaches x_wec

    Nx = 501

    #
    # initialization
    #
    def __init__(self,
            surfDir=surfDir,
            MLERinput=MLERinput,
            dt=dt,
            toffset=toffset,
            x_wec=x_wec,
            x_inlet=x_inlet):

        # inputs
        self.surfDir    = surfDir
        self.MLERinput  = MLERinput
        self.dt         = dt
        self.toffset    = toffset
        self.x_wec      = x_wec
        self.x_inlet    = x_inlet

        # outputs
        self.csvfiles = []
        self.indices = []
        self.ipeak = -1

        self.times  = []
        self.z_sim  = []
        self.z_ref  = []
        self.z0_sim = []
        self.z0_ref = []

        self.tpeak = -1
        self.xpeak = []
        self.zpeak = []
        self.xpeak_ref  = []
        self.zpeak_ref  = []

    #
    # find csv files
    #
    def find_csv(self,surfDir=''):
        if surfDir=='': surfDir = self.surfDir
        self.csvfiles = []
        for f in os.listdir(surfDir):
            if f.endswith(".csv"):
                self.csvfiles.append(surfDir+os.sep+f)
        self.Ntimes = len(self.csvfiles)
        return self.csvfiles

    #
    # sort csv files
    #
    def sort_csv(self):
        self.times = np.zeros((self.Ntimes))
        i = -1
        nearpeak = 9e9
        for csv in self.csvfiles:
            #t = float( csv[:-4].split('_')[-1] )
            t = float( csv[:-4].split('_')[-1] ) * self.dt
            i += 1
            self.times[i] = t
            delta = np.abs(t-self.toffset)
            if delta < nearpeak:
                nearpeak = delta
                self.ipeak = i
        self.indices = [ i[0] for i in sorted(enumerate(self.times), key=lambda x:x[1]) ]
        self.times   = self.times[self.indices]
        return self.times, self.indices

    #
    # define csv reader
    #
    def read_csv(self,fname,N=-1): # assume 1 header line
        if N < 0: 
            with open(fname,'r') as f:
                for line in f: N+=1
        x = np.zeros((N))
        y = np.zeros((N))
        with open(fname,'r') as csv:
            i = -1
            csv.readline()
            for line in csv:
                i += 1
                line = line.split(',')
                x[i] = float(line[0])
                y[i] = float(line[1])
        order = x.argsort()
        x = x[order]
        y = y[order]

        # from post_wave_fft.py
        #if any(np.diff(x)==0):
        if any(np.abs(np.diff(x))<1e-6):
            print 'Warning: duplicate x found in',fname
# SKIP FIX# {{{
#            newx = np.zeros((N))
#            newy = np.zeros((N))
#            newN = N
#            i = 0
#            inew = 0
#            while i < N:
#                if i < N-1 and x[i+1]==x[i]:
#                    i2 = i+2
#                    while i2 < N and x[i2]==x[i]:
#                        i2 += 1
#                    newx[inew] = x[i]
#                    newy[inew] = np.mean(y[i:i2])
#                    ndup = i2-i-1
#                    i += ndup+1
#                    newN -= ndup
#                else: 
#                    newx[inew] = x[i]
#                    newy[inew] = y[i]
#                    i+= 1
#                inew += 1
#            x = newx[:newN]
#            y = newy[:newN]
# }}}
            newx = np.zeros((N))
            newy = np.zeros((N))
            norm = np.ones((N))
            newN = N
            inew = 0
            newx[0] = x[0]
            newy[0] = y[0]
            for i in range(1,N):
                if x[i] - x[inew] < 1e-6: #duplicate x
                    print 'DUPLICATE x ({:f},{:f}) == ({:f},{:f})'.format(x[inew],y[inew],x[i],y[i])
                    newN -= 1
                    newy[inew] += newy[i]
                    norm[inew] += 1
                else:
                    inew += 1
                    newx[inew] = x[i]
                    newy[inew] = y[i]
            print '  number of duplicated points :',np.count_nonzero(norm>1)
            print '  maximum number overlapping  :',np.max(norm)
            x = newx[:newN]
            y = newy[:newN] / norm[:newN]
            print '  old length',N,'new length',newN,len(x)

        return x,y

    #
    # process all csv files
    #
    def process_all_times(self):
        self.z_sim   = np.zeros((self.Ntimes))
        self.z0_sim  = np.zeros((self.Ntimes))
        read_data = True

        try:
            tread = np.zeros((self.Ntimes))
            with open('trace.dat','r') as f:
                print 'Attempting read of existing trace.dat'
                i = 0
                f.readline() #header
                for line in f:
                    vals = [ float(val) for val in line.split() ]
                    tread[i] = vals[0]
                    self.z_sim[i] = vals[1]
                    self.z0_sim[i] = vals[2]
                    i += 1
            with open('peak.dat','r') as f:
                print 'Attempting read of existing peak.dat'
                f.readline() #header
                x,z = [],[]
                for line in f:
                    vals = [ float(val) for val in line.split() ]
                    x.append( vals[0] )
                    z.append( vals[1] )
                self.xpeak = np.array(x)
                self.zpeak = np.array(z)

            if np.max(np.abs( tread-self.times )) < 1e-12:
                self.tpeak = self.times[self.ipeak]
                read_data = False
            else:
                print 'Different times -- read:',tread
                print '  expected:',self.times
                print '  max diff:',np.max(np.abs( tread-self.times ))

        except IOError: pass

        if read_data:
            for itime in range(self.Ntimes):
                idx     = self.indices[itime]
                t       = self.times[itime]
                fname   = self.csvfiles[idx]

                print 't= {:f} s : {:s}'.format(t,fname)
                x,z = self.read_csv(fname)

                self.z_sim[itime]  = np.interp(self.x_wec,x,z)
                self.z0_sim[itime] = np.interp(self.x_inlet,x,z)

                if itime==self.ipeak:
                    self.xpeak = x
                    self.zpeak = z
                    self.tpeak = t
            # end of time loop

            # dump out trace at device location
            with open('trace.dat','w') as f:
                f.write('#Data: t z(x=0,t) z(x=x_in,t)\n')
                for ti,zi,z0i in zip(self.times,self.z_sim,self.z0_sim):
                    f.write('%f %f %f\n'%(ti,zi,z0i))

            # dump out peak profile
            with open('peak.dat','w') as f:
                f.write('#Data: x z(x,t=0)\n')
                for xi,zi in zip(self.xpeak,self.zpeak):
                    f.write('%f %f\n'%(xi,zi))

        # endif read_data

        print 'Max z(x=0,t)=',np.max(self.z_sim),'at t=',self.times[np.argmax(self.z_sim)]

        # numerically evaluate the derivative

#        dx = self.xpeak[1:] - self.xpeak[:-1]
#        dzdx = (self.zpeak[1:] - self.zpeak[:-1]) / dx
#        dzdx = (self.zpeak[1:] - self.zpeak[:-1]) / dx

#        dx = np.diff( self.xpeak )
#        if any( dx==0 ): print 'Warning diff(xpeak)==0'
#        for i in range(1,len(self.xpeak)):
#            if np.abs(self.xpeak[i] - self.xpeak[i-1]) < 1e-6:
#                print 'DUPLICATE ({:f},{:f}) ({:f},{:f})'.format( 
#                        self.xpeak[i-1], self.zpeak[i-1], self.xpeak[i], self.zpeak[i] )
#        dz = np.diff( self.zpeak )
#        dzdx = dz / dx
#        imax = np.argmax(dzdx)
#        x_dzdx_max = 0.5*(self.xpeak[imax]+self.xpeak[imax+1])
#        dzdx_max = dzdx[imax]
        
        dx = self.xpeak[2:] - self.xpeak[:-2]
        dz = self.zpeak[2:] - self.zpeak[:-2]
        dzdx = dz / dx
        imax = np.argmax(np.abs(dzdx))
        x_dzdx_max = self.xpeak[imax+1]
        dzdx_max = dzdx[imax]

        print 'Approx max dz(x,tpeak)/dx=',dzdx_max,'at x=',x_dzdx_max

        return self.xpeak, self.zpeak, self.z_sim, self.z0_sim

    #
    # calculate reference values
    #
    def calc_ref(self):
        w = [] # frequency, rad/s
        S = [] # spectral amplitude, m^2
        p = [] # phase, rad
        k = [] # wavenumber, rad^2/m
        print 'Reading MLER input file:',self.MLERinput
        print '--- begin header ---'
        with open(self.MLERinput,'r') as f:
            for line in f:
                if line.strip()=='': break
                if line.startswith('#'):
                    sys.stdout.write(line)
                    continue
                wi,Si,pp,ki = [ float(val) for val in line.split() ]
                w.append(wi)
                S.append(Si)
                p.append(pp)
                k.append(ki)
        print '--- end header ---'

        self.w  = np.array(w)
        self.dw = w[1] - w[0] #dw = np.diff(w)
        self.S  = np.array(S)
        self.A  = np.array(2*self.S*self.dw)**0.5 # updated 6/7/16 after spectral amplitude was corrected in final version
        self.p  = -np.array(p) # updated 2/10/16 after *coeffs*.txt output was changed
        self.k  = np.array(k)

        self.z_ref  = np.zeros((self.Ntimes))
        self.z0_ref = np.zeros((self.Ntimes))
        t = self.times - self.toffset
        for i in range(self.Ntimes):
            # updated 6/7/16 after spectral amplitude was corrected in final version
            self.z_ref[i]  = np.sum( self.A*np.cos( self.k*self.x_wec   - self.w*t[i] + self.p ) )#*self.dw
            self.z0_ref[i] = np.sum( self.A*np.cos( self.k*self.x_inlet - self.w*t[i] + self.p ) )#*self.dw

        print 'Calculating wave profile',self.ipeak,'at t =',self.tpeak,'~= tpeak'
        self.xpeak_ref  = np.linspace(-600,600,self.Nx)
        xref = self.xpeak_ref - self.x_wec
        self.zpeak_ref  = np.zeros((self.Nx))
        for i in range(self.Nx):
            # updated 6/7/16 after spectral amplitude was corrected in final version
            self.zpeak_ref[i] = np.sum( self.A*np.cos( self.k*xref[i] - self.w*t[self.ipeak] + self.p ) )#*self.dw

        # numerically evaluate the derivative
##        dx = self.xpeak_ref[1:] - self.xpeak_ref[:-1]
##        dzdx = (self.zpeak_ref[1:] - self.zpeak_ref[:-1]) / dx
#        dx = np.diff( self.xpeak_ref )
#        if any( dx==0 ): print 'Warning diff(xpeak_ref)==0'
#        dz = np.diff( self.zpeak_ref )
#        dzdx = dz / dx
#        imax = np.argmax(dzdx)
#        x_dzdx_max = 0.5*(self.xpeak_ref[imax]+self.xpeak_ref[imax+1])
#        dzdx_max = dzdx[imax]
#        print 'Approx max ref dz(x,tpeak)/dx=',dzdx_max,'at x=',x_dzdx_max

        # analytically evaluate the derivative at t=0 (peak)
        dzdx_ref  = np.zeros((len(self.xpeak)))
        for i,x in enumerate(self.xpeak):
            dzdx_ref[i] = np.sum( -self.A*self.k*np.sin( self.k*x + self.p ) )*self.dw
        imax = np.argmax(np.abs(dzdx_ref))
        x_dzdx_max = self.xpeak[imax]
        dzdx_max = dzdx_ref[imax]

        print 'Max ref dz(x,tpeak)/dx=',dzdx_max,'at x=',x_dzdx_max

        return self.z_ref, self.z0_ref, self.xpeak_ref, self.zpeak_ref

    #
    # comparison plots
    #
    def compare_hist(self):
        fig, (ax0,ax1) = plt.subplots(2,sharex=True)
        plt.subplots_adjust(hspace=0.25) #default hspace=0.2

        ax0.plot(self.times,self.z0_ref,'k-',label='MLER input')
        ax0.plot(self.times,self.z0_sim,'b.-',label='Star sim')
        ax0.set_ylabel('z( x={:.2f}, t )'.format(self.x_inlet))
        ax0.set_title('Wave height at inlet location')
        ax0.legend(loc='best')

        ax1.plot(self.times,self.z_ref,'k-',label='MLER input')
        ax1.plot(self.times,self.z_sim,'b.-',label='Star sim')
        ax1.set_xlabel('t')
        ax1.set_ylabel('z( x= {:.2f}, t )'.format(self.x_wec))
        ax1.set_title('Wave height at device location')

        fig.savefig('focusedwave_trace.png')

    def compare_peak_profile(self):
        fig, ax = plt.subplots(1)
        ax.plot(self.xpeak_ref,self.zpeak_ref,'k-',label='MLER input')
        ax.plot(self.xpeak,self.zpeak,'b.-',label='Star sim')
        ax.set_xlabel('x')
        ax.set_ylabel('z( x, tpeak )')
        ax.set_title('Focused wave peak (t = {:f} s)'.format(self.tpeak))
        ax.legend(loc='best')

        fig.savefig('focusedwave_peak_profile.png')

#===================================================================
if __name__ == '__main__':

    if len(sys.argv) > 1: surfDir = sys.argv[1]
    fw = focused_wave()

    csvfiles        = fw.find_csv()
    times, indices  = fw.sort_csv()
    xpeak, zpeak, z_sim, z0_sim = fw.process_all_times()

    z_ref, z0_ref, xpeak_ref, zpeak_ref = fw.calc_ref()

    fw.compare_hist()
    fw.compare_peak_profile()
    plt.show()

