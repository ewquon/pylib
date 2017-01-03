#!/usr/bin/python
"""
Look for a file called "refs.yaml" and load reference values to be used during
any post processing.
"""
import os
import yaml

def read( ref=dict(), refFile='refs.yaml', verbose=False ):
    """ Read reference values, updating an optional input dictionary with default values
    """
    if os.path.isfile(refFile):
        refDict = yaml.load(open(refFile,'r'))
        for key,val in refDict.iteritems():
            if verbose:
                if key in ref:
                    print 'Overwriting',key,ref[key],'with',val
                else:
                    print 'Setting',key,'to',val
            ref[key] = val
        if verbose:
            print 'Reference values:',ref
    else:
        if verbose: print refFile,'containing reference values not found'
    return refData(ref)

class refData(object):
    def __init__(self,dataDict):
        self._dict = dataDict
        for key,val in dataDict.iteritems():
            setattr(self,key,val)

    def __str__(self):
        s = ''
        for key,val in self._dict.iteritems():
            s += '\n  {:s} : {:}'.format(key,val)
        return s

    def __bool__(self):
        return len(self._dict) > 0

    def __nonzero__(self):
        return self.__bool__()

    def dict(self):
        return self._dict

#------------------------------------------------------------------------------
if __name__ == '__main__':
    refValues = read({'Uref':10.0},verbose=True)
