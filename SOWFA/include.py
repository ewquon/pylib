#!/usr/bin/env python
#
# Scrape parameters from specified SOWFA include file
#
from __future__ import print_function
import os

def readVarDefinition(line):
    if '//' in line:
        line = line.split('//')
        comment = line[1].strip()
        part1 = line[0].strip()
    else:
        comment = ''
        part1 = line.strip()
    assert( part1[-1]==';' )

    nameval = part1[:-1].split()
    name = nameval[0]
    val = ' '.join(nameval[1:])
    try: 
        val = float(val)
    except ValueError:
        pass

    return (name,val,comment)

def read(fname='setUp',verbose=False):
    params = dict()
    comments = dict()
    if not os.path.isfile(fname):
        print('File not found:',fname)
        return params

    with open(fname,'r') as f:
        processLine = True
        for line in f:
            if line.startswith('//') or line.startswith('#') or line.strip()=='':
                continue
            elif line.startswith('/*'):
                processLine = False
            elif not processLine and '*/' in line:
                processLine = True
            elif processLine:
                p = readVarDefinition(line)
                if not p is None:
                    if verbose: print(p)
                    params[p[0]] = p[1]
                    comments[p[0]] = p[2]

    return params, comments

def read_setUp(fpath=None):
    if fpath is None:
        fpath = 'setUp'
        # find path to setUp file
        while not os.path.isfile(fpath):
            fpath = '..' + os.sep + fpath
            if not os.path.isdir(os.path.dirname(fpath)):
                break
    setUp = dict()
    if not os.path.isfile(fpath):
        print('setUp file not found')
    else:
        with open(fpath,'r') as f:
            for _ in range(7):
                f.readline() # skip foam header
            for line in f:
                line = line.strip()
                if line.startswith('//'):
                    continue
                line = line.split()
                if len(line) <= 1:
                    continue
                if line[1].endswith(';'):
                    val = line[1][:-1]
                    try:
                        val = float(val)
                    except ValueError:
                        pass
                    else:
                        if val.is_integer(): val = int(val)
                    setUp[line[0]] = val
    return setUp
