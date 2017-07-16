#
# Helper module for scripts producing lots of repetitive screen output.
# 
# Written by Eliot Quon (eliot.quon@nrel.gov) -- 2017-07-16
#
import sys

class screen(object):

    def __init__(self,initString=''):
        if len(initString) > 0:
            sys.stdout.write(initString)
        self.last = initString
        self.N = len(self.last)

    def overwrite(self,s):
        N = len(s)
        if N < self.N:
            newstr = s + (self.N-N)*' '
        else:
            newstr = s
        sys.stdout.write('\r'+newstr)
        sys.last = s
        self.N = N
    # alias:
    out = overwrite

    def endl(self):
        sys.stdout.write('\n')
        self.N = 0
