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
    dataobj = refData(ref)
    dataobj.read(refFile)
    return dataobj

class refData(object):
    def __init__(self,dataDict):
        self._dict = dataDict
        self._updateData()

    def _updateData(self):
        for key,val in self._dict.iteritems():
            setattr(self,key,val)

    def __getattribute__(self,name):
        try:
            return object.__getattribute__(self,name)
        except AttributeError:
            return 'n/a'

    def __str__(self):
        s = ''
        for key,val in self._dict.iteritems():
            s += '\n  {:s} : {:}'.format(key,val)
        return s

    def __bool__(self):
        return len(self._dict) > 0

    def __nonzero__(self):
        return self.__bool__()

    def dict(self): return self._dict

    def read(self,refFile,verbose=False):
        if not os.path.isfile(refFile):
            if verbose: print refFile,'containing reference values not found'
            return

        refDict = yaml.load(open(refFile,'r'))
        if refDict is None:
            if verbose: print 'No reference values read from',refFile
            return

        for key,val in refDict.iteritems():
            if verbose:
                if key in self._dict:
                    print 'Overwriting',key,self._dict[key],'with',val
                else:
                    print 'Setting',key,'to',val
            self._dict[key] = val

        if verbose: print 'Reference values:',ref
        self._updateData()


    def readABLProperties(self,fname='constant/ABLProperties',verbose=False):
        with open(fname,'r') as f:
            skip = False
            parensLevel = 0
            bracesLevel = 0
            for line in f:
                if not skip and '/*' in line:
                    skip = True
                elif skip and '*/' in line:
                    skip = False
                elif not skip and ( '(' in line or ')' in line or '{' in line or '}' in line ):
                    # only read single value lines, not vectors/arrays
                    parensLevel += line.count('(')
                    parensLevel -= line.count(')')
                    # don't read dictionaries either...
                    bracesLevel += line.count('{')
                    bracesLevel -= line.count('}')
                elif not skip and parensLevel==0 and bracesLevel==0:
                    line = line.strip()
                    if not line.startswith('//') and ';' in line:
                        line = line.split(';')
                        comment = line[1].strip()
                        parts = line[0].split()
                        var = parts[0]
                        val = ' '.join(parts[1:])
                        val = val.replace('"','').strip()
                        if not val=='' and not val[0]=='$':
                            if verbose:
                                if var in self._dict: print 'ABLProperties: Overwriting',var,self._dict[var],'with',val
                                else: print 'ABLProperties: Setting',var,'to',val
                            self._dict[var] = val
            # done looping over lines
        # finished reading file
        self._updateData()


#------------------------------------------------------------------------------
if __name__ == '__main__':
    refValues = read({'Uref':10.0},verbose=True)
