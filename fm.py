#!/usr/local/bin/python

class fm:
    def __init__(self,Fp,Fv):
        self.p = []
        self.v = []
        self.total = []
        self.addTime(Fp,Fv)

    def addTime(self,Fp,Fv):
        #for Fpi,Fvi in zip(Fp,Fv):
        #    print Fpi,Fvi
        #    self.p.append(Fpi)
        #    self.v.append(Fvi)
        #    self.total.append(Fpi+Fvi)
        self.p.append(Fp)
        self.v.append(Fv)
        self.total.append(Fp+Fv)
