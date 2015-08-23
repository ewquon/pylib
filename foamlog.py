#!/usr/bin/python

class foamLog:
    """Class to parse OpenFOAM log files"""

    def __init__( self, fname,\
            timeStr='Time = ',\
            solveStr='Solving for ',\
            residStr='Final residual = ',\
            #forceStrs=['Cl', 'Cd', 'Cm' ] ):
            #forceStrs=['Cl    =', \
            #           'Cd    =', \
            #           'Cm    =' ] ):
            aeroStrs=['Cl', 'Cd' ] ):
        self.fname = fname
        self.timeStr = timeStr
        self.solveStr = solveStr
        self.residStr = residStr
        self.aeroStrs = aeroStrs

        self.it = []
        self.t = []
        self.R = dict()
        self.C = dict()
        self.F = dict()
        self.varlist = []

# TODO: CLEAN UP THIS IMPLEMENTATION, AND/OR READ FROM postProcessing/*
    def read(self):
        N=0
        readForces = False
        readMoments = False
        with open(self.fname,'r') as f:
            for line in f:

                if line[:len(self.timeStr)]==self.timeStr: 
                    self.t.append(float(line[len(self.timeStr):]))
                    N += 1
                    self.it.append(N)

                elif self.solveStr in line:
                    varline = line.split(self.solveStr)[1]
                    var = varline.split(',')[0]

                    valline = line.split(self.residStr)[1]
                    val = float(valline.split(',')[0])

                    try: 
                        if len(self.R[var]) < N:
                            self.R[var].append(val)
                        elif len(self.R[var])==N:
                            self.R[var][N-1] = val
                        else:
                            sys.exit('too many values in residual list')
                    except KeyError: 
                        self.varlist.append(var)
                        self.R[var] = [val]

                elif line.strip().startswith('sum of forces:'): 
                    if len(self.F)==0:
                        self.F['Fx'] = []
                        self.F['Fy'] = []
                        self.F['Fz'] = []
                        self.F['Mx'] = []
                        self.F['My'] = []
                        self.F['Mz'] = []
                    readForces = True
                elif line.strip().startswith('sum of moments:'): 
                    readMoments = True
                elif readForces or readMoments:
                    line = line.strip()
                    if line.startswith('pressure'):
                        line = line.split('(')[1][:-1].split()
                        if readForces:
                            [Fx,Fy,Fz] = [float(s) for s in line]
                        elif readMoments:
                            [Mx,My,Mz] = [float(s) for s in line]
                    elif line.startswith('viscous'):
                        line = line.split('(')[1][:-1].split()
                        if readForces:
                            self.F['Fx'].append( Fx + float(line[0]) )
                            self.F['Fy'].append( Fy + float(line[1]) )
                            self.F['Fz'].append( Fz + float(line[2]) )
                            readForces = False
                        elif readMoments:
                            self.F['Mx'].append( Mx + float(line[0]) )
                            self.F['My'].append( My + float(line[1]) )
                            self.F['Mz'].append( Mz + float(line[2]) )
                            readMoments = False

                elif any([line.strip().startswith(substr) for substr in self.aeroStrs]):

                    line = line.split()
                    var = line[0]
                    val = float(line[-1])

                    try:
                        if len(self.C[var]) < N:
                            self.C[var].append(val)
                        else:
                            sys.exit('too many values in q array')
                    except KeyError:
                        self.C[var] = [val]

        for rv in self.R: # reduce N if we have an incomplete log file
            if len(self.R[rv]) < N:
                print len(self.R[rv]),'<',N,': trimming N for incomplete log file'
                N -= 1
                break
        for rv in self.F: # reduce N if we have an incomplete log file
            if len(self.F[rv]) < N:
                print len(self.F[rv]),'<',N,': trimming N for incomplete log file'
                N -= 1
                break
        for rv in self.C: # reduce N if we have an incomplete log file
            if len(self.C[rv]) < N:
                print len(self.C[rv]),'<',N,': trimming N for incomplete log file'
                N -= 1
                break
        for rv in self.R:
            self.R[rv] = self.R[rv][:N]
        for fv in self.F:
            self.F[fv] = self.F[fv][:N]
        for fv in self.C:
            self.C[fv] = self.C[fv][:N]
        print 'Read',N,'steps:'
        print '  residuals  ',self.R.keys()
        print '  forces     ',self.F.keys()
        print '  aero coeffs',self.C.keys()

        return self.t[:N], self.R, self.F, self.C

    def Rvars(self):
        return self.R.keys()

    def Fvars(self):
        return self.F.keys()

    def Cvars(self):
        return self.C.keys()

