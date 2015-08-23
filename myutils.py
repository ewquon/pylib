#!/usr/bin/python
import numpy as np

#-------------------------------------------------------------------------------
def smartTime(t):
    threshold = 100.0;
    tUnit = ['s','min','hrs','days']
    tConv = [60,60,24]
    i = 0
    t *= 1.0 # convert to float
    while t > threshold and i < len(tUnit)-1:
        t /= tConv[i]
        i += 1
    return '{:.1f} {}'.format(t,tUnit[i])

#-------------------------------------------------------------------------------
def read_generic_data(fname,delim=None,headerLines=0,cast=None):
#    data = []
#    with open(fname,'r') as f:
#        for i in range(headerLines): f.readline()
#        line = f.readline().split(delim)
#        Nvars = len(line)
#        for val in line: data.append([val])
#        for line in f:
#            line = line.split(delim)
#            for i in range(Nvars):
#                data[i].append(line[i])
#    if not cast:
#        return data
#    else:
#        formattedData = []
#        for col in data:
#            formattedData.append([cast(val) for val in col])
#        return formattedData
    with open(fname,'r') as f:
        for i in range(headerLines): f.readline()
        line = f.readline().split(delim)
        Nvars = len(line)
        Ndata = 0
        for line in f: Ndata += 1
    data = np.zeros((Ndata,Nvars))
    with open(fname,'r') as f:
        for i in range(headerLines): f.readline()
        for i in range(Ndata):
            rowvals = f.readline().split(delim)
            if cast:
                data[i,:] = [ cast(val) for val in rowvals ]
            else:
                data[i,:] = rowvals
    return data

#-------------------------------------------------------------------------------
def floatList(L): return [float(s) for s in L]

#-------------------------------------------------------------------------------
# DOES THIS WORK?
def parse_openfoam(s):
    if not '(' in s and not ')' in s: return s.split()
    slist = s.split()
    newList = []
    curx = ''
    for x in slist:
        curx += ' ' + x
        #print 'curx:',curx
        if curx.count('(')==curx.count(')'):
            curx = curx.strip()
            if curx.startswith('('): curx = curx[1:-1]
            #print 'attempt to add', curx
            newList.append( parse_openfoam(curx) ) # RECURSE
            curx = ''
            #print 'current list:',newList
    return newList

#-------------------------------------------------------------------------------
def moving_avg(u,window_size):
    newu = u
    N = len(u)
    for i in range(N):
        ist = i-window_size
        ind = i+window_size
        if ist < 0:
            ind -= ist
            ist = 0
        elif ind >= N:
            ist -= (ind-N+1)
            ind = N-1
        newu[i] = np.sum(u[ist:ind+1])/(2*window_size+1)
    return newu

def moving_median(u,window_size):
    newu = u
    N = len(u)
    for i in range(N):
        ist = i-window_size
        ind = i+window_size
        if ist < 0:
            ind -= ist
            ist = 0
        elif ind >= N:
            ist -= (ind-N+1)
            ind = N-1
        newu[i] = np.median(u[ist:ind+1])
    return newu

def moving_op(fn,dataset,window_size):
    dataset = list(dataset)
    for u in enumerate(dataset):
        #newu = u[1] # newu is a reference to u[1]!!!
        newu = np.zeros(u[1].shape)
        N = len(newu)
        for i in range(N):
            ist = i-window_size
            ind = i+window_size
            if ist < 0:
                ind -= ist
                ist = 0
            elif ind >= N:
                ist -= (ind-N+1)
                ind = N-1
            newu[i] = fn(u[1][ist:ind+1])
        dataset[u[0]] = newu
    return dataset


#-------------------------------------------------------------------------------

def isfloat(x):
    try:
        a = float(x)
    except ValueError:
        return False
    else:
        return True

def isint(x):
    try:
        int(x)
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b


