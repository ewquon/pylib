#!/usr/local/bin/python
# Works with blade_def class
import numpy as np

class BEMT:
    TOL = 1e-6
    MAXITER = 100
    VERBOSE = True

    def __init__(self,blade,dr,R=1.0,B=3,rhoinf=1.0,Uinf=1.0,omega=1.0,theta0=0.0,tipCorrection=''):
        self.B = B                          # number of blades
        self.R = R                          # blade radius, m
        self.rhoinf = rhoinf                # freestream velocity, m/s
        self.Uinf = Uinf                    # freestream velocity, m/s
        self.omega = omega*np.pi/30         # rotor rotation rate, RPM
        self.theta0 = theta0*np.pi/180      # collective pitch, deg (note: convention is +ve toward feather)
        self.blade = blade                  # blade_def object containing geometry
        self.dr = dr                        # length of each strip
        self.tipCorrection = tipCorrection  # prandtl, xusankar

        self.A = np.pi * R**2
        self.r = np.array([ section.r for section in blade.sections ])
        self.N = len(blade.sections)


    def calculate(self):

        self.Cl = np.zeros(self.N)
        self.Cd = np.zeros(self.N)
        self.F = np.ones(self.N)

        self.sig = self.B/(2*np.pi*self.r) * np.array([ section.chord for section in self.blade.sections ])

        #self.Utip = self.omega*self.R
        #self.TSR = self.Utip / self.Uinf
        #if self.VERBOSE: print 'Tip speed ratio :',self.TSR

        local_inflow_ratio = self.omega * self.r / self.Uinf

        # note: convention is +ve toward feather
        self.theta = self.theta0 + np.pi/180.0*np.array([ section.twist for section in self.blade.sections ])
        #if self.VERBOSE: print 'Local pitch :',self.theta*180/np.pi

        #
        # Solve for induction factors with iterative approach (Manwell p116)
        #

        # 1. initial guesses for a and a'
        self.a = np.ones(self.N)/3 # Betz limit
        self.aswirl = np.zeros(self.N) #no swirl

        change1 = self.TOL + 1
        change2 = self.TOL + 1
        niter = 0
        while change1 > self.TOL and change2 > self.TOL and niter < self.MAXITER:
            niter += 1
            last_a = self.a
            last_aswirl = self.aswirl

            # 2. calculate angle of relative wind
            self.phi = np.arctan( (1-self.a) / ((1+self.aswirl)*local_inflow_ratio) ) # local inflow angle
            if self.tipCorrection: self.update_tip_correction()

            # 3. update aero
            self.alpha = self.phi - self.theta
            for i in range(self.N):
                self.Cl[i], self.Cd[i], Cm = self.blade.sections[i].aero_lookup( self.alpha[i] )

            # 4. update a, a'
            self.a = 1 / ( 1 + 4*self.F*np.sin(self.phi)**2/(self.sig*self.Cl*np.cos(self.phi) ) )
            self.aswirl = self.a * np.tan(self.phi) / local_inflow_ratio
            #print 'a ',self.a
            #print 'a\'',self.aswirl
            change1 = np.max(np.abs(self.a - last_a))
            change2 = np.max(np.abs(self.aswirl - last_aswirl))
            if self.VERBOSE: print ' *** iter',niter,': max change in a, a\' =',change1,change2,'***'

        if self.VERBOSE:
            #if self.tipCorrection:
            #    print self.tipCorrection,'tip loss correction factor :',self.F
            #print 'a converged to',self.a
            #print 'a\' converged to',self.aswirl
            print 'max induction factor :',np.max(self.a)
        assert(np.max(self.a) <= 0.5)

        #norm = self.sig * np.pi * self.rhoinf * (self.Uinf*(1-self.a))**2 / np.sin(self.phi)**2
        #self.dFN = norm*( self.Cl*np.cos(self.phi) + self.Cd*np.sin(self.phi) )*self.r
        #self.dQ  = norm*( self.Cl*np.sin(self.phi) - self.Cd*np.cos(self.phi) )*self.r**2
        #self.T = np.dot( self.dFN, self.dr )
        #self.Q = np.dot( self.dQ,  self.dr )
        #print 'Integrated thrust [kN], torque [MW]:', self.T/1000, self.Q/1e6

        dT = self.F * self.rhoinf * self.Uinf**2 * 4*self.a*(1-self.a) * np.pi * self.r
        dQ = 4*self.F*self.aswirl*(1-self.a) * self.rhoinf*self.Uinf * np.pi * self.r**3 * self.omega
        print 'Integrated thrust [kN], power [MW]:', np.dot(dT,self.dr)/1000, self.omega*np.dot(dQ,self.dr)/1e6


    def update_tip_correction(self):
        #self.F[:] = 1.0
        #if self.tipCorrection.lower() == 'prandtl':
        f = self.B/2*( (self.R-self.r)/(self.r*self.phi) )
        self.F = 2/np.pi * np.arccos(np.exp(-f))

        # correction from Xu and Sankar 2002
        if self.tipCorrection.lower() == 'xusankar':
            threshold = 0.7*self.R
            for i in range(self.N):
                if self.r[i] < threshold:
                    self.F[i] = 1 - (1-self.F[i])*self.r[i]/threshold
                else:
                    self.F[i] = (self.F[i]**0.85 + 0.5) / 2

    def print_results_table(self):
        row_format = "{:>8}" + "{:>17}" * 7
        print row_format.format('station','blade_pitch','a','a\'','F','angle_of_attack','Cl','Cd')
        for i in range(self.N):
            print row_format.format(i+1, self.theta[i]*180.0/np.pi, self.a[i], self.aswirl[i], self.F[i], \
                    self.alpha[i]*180.0/np.pi, self.Cl[i], self.Cd[i])

