#
# RBF interpolation module
# written by Eliot Quon (eliot.quon@nrel.gov)
# 
from __future__ import print_function
import numpy as np

def tps(r):
    """Thin-plate spline RBF"""
    return np.array([ ri*ri*np.log(ri) if ri > 0 else 0.0 for ri in r ])

class RBFInterpolant(object):
    """An instance of the interpolant is associated with a given set of
    data sites"""

    def __init__(self, x, fvals=None,
                 func=tps, P=lambda x: [1, x[0], x[1]]):
        """Data sites (or "centers") with shape (N,dim)
        Default RBF is the thin-plate spline
        Default polynomial is P(x,y) = B0 + B1*x + B2*y
        """
        self.x = x
        self.N, self.dim = x.shape
        self.P = P
        self.Q = len(P(x[0,:]))
        self.func = func
        self.calculate_LHS()
        if fvals is not None:
            self.update(fvals)
        else:
            self.Ntimes = 0
            self.coefs = None

    def calculate_LHS(self):
        """Populate the left-hand side of the interpolation system
        for a given set of data sites and specified RBF"""
        Amat = np.zeros((self.N+self.Q, self.N+self.Q))
        for j in range(self.N):
            r2 = np.zeros((self.N,))
            for d in range(self.dim):
                r2 += (self.x[j,d] - self.x[:,d])**2
            Amat[j,:self.N] = self.func(np.sqrt(r2))
            pvals = self.P(self.x[j,:])
            Amat[j,self.N:] = pvals
            Amat[self.N:,j] = pvals
        self.A = Amat

    def update(self,fvals):
        """Calculates interpolation coefficients for a set of function
        values associated with the data sites. Values should have shape
        (N,Ntimes)"""
        N,self.Ntimes = fvals.shape
        assert(N == self.N)
        coefs = np.zeros((self.N+self.Q,self.Ntimes))
        RHS = np.zeros(np.shape(coefs))
        RHS[:N,:] = fvals
        for itime in range(self.Ntimes):
            coefs[:,itime] = np.linalg.solve(self.A, RHS[:,itime])
        self.coefs = coefs

    def evaluate(self,xi):
        """Evaluates coefficients at xo for all times"""
        assert(self.coefs is not None)
        r2 = np.zeros((self.N,))
        interpvec = np.zeros((self.N+self.Q,))
        for d in range(self.dim):
            r2 += (xi[d] - self.x[:,d])**2
        interpvec[:self.N] = self.func(np.sqrt(r2))
        interpvec[self.N:] = self.P(xi)
        return np.dot(interpvec, self.coefs)
