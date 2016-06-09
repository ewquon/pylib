#!/usr/bin/python
import sys
import os
import numpy as np

class csv:

    def __init__(self,fname='',outputnames=[]):
        self.fname = fname
        self.outputnames = outputnames
        self.Noutputs = -1
        self.x = []
        self.y = []

    def find_latest(self,prefix='',searchDir='.'):
        csvfiles = []
        for f in os.listdir(searchDir):
            if f.startswith(prefix) and f.endswith(".csv"):
                csvfiles.append(searchDir+os.sep+f)
        times = np.zeros((len(csvfiles)))
        i = -1
        for csv in csvfiles:
            i += 1
            times[i] = float( csv[:-4].split('_')[-1] )
        idx = np.argmax(times)
        self.fname = csvfiles[idx]
        return self.fname

    def read(self,fname='',N=-1,sort=False):
        if fname=='': fname = self.fname
        if N < 0:
            N = 0
            with open(fname,'r') as f:
                hdr = f.readline().split(',')
                for line in f: N+=1
        Noutputs = len(hdr) - 1 # don't count the independent variable
        self.Noutputs = Noutputs
        if len(self.outputnames)==0:
            self.outputnames = [ s.replace('"','') for s in hdr ]
        else: assert( len(self.outputnames)==Noutputs )
        x = np.zeros((N))
        y = np.zeros((N,Noutputs))
        with open(fname,'r') as f:
            i = -1
            f.readline()
            for line in f:
                i += 1
                line = line.split(',')
                x[i] = float(line[0])
                y[i,:] = [ float(val) for val in line[1:] ]
        self.x = x
        self.y = y
        if sort: self.sort()
        return self.x, self.y

    def sort(self):
        order = self.x.argsort()
        self.x = self.x[order]
        for j in range(Noutputs):
            self.y[:,j] = self.y[order,j]

