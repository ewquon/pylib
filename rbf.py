#
# RBF interpolation module
# written by Eliot Quon (eliot.quon@nrel.gov)
# 
from __future__ import print_function
import sys
import numpy as np

def thin_plate(r):
    """Thin-plate spline RBF"""
    return np.array([ ri*ri*np.log(ri) if ri > 0 else 0.0 for ri in r ])
def grad_thin_plate(r):
    return np.array([ 2*np.log(ri) + 1.0 if ri > 0 else 0.0 for ri in r ])

def multiquadric(r,eps=1.0):
    """Multiquadric RBF, c > 0"""
    return np.sqrt(1 + (eps*r)**2)
def grad_multiquadric(r,eps=1.0):
    return eps**2*r / np.sqrt(1 + (eps*r)**2)

def inverse(r,eps=1.0):
    """Inverse quadratic RBF, c > 0"""
    return 1.0 / (1 + (eps*r)**2)
def grad_inverse(r,eps=1.0):
    return -2*eps**2*r / (1 + (eps*r)**2)**2

def gaussian(r,eps=1.0):
    """Gaussian RBF"""
    return np.exp(-(eps*r)**2)
def grad_gaussian(r,eps=1.0):
    return -2*eps**2*r * np.exp(-(eps*r)**2)


class RBFInterpolant(object):
    """Radial basis function interpolant. 
    
    An instance of the interpolant is associated with a given set of
    data sites (or "centers"). The workflow for interpolating a set of
    data (x,f(x)) to point xi is as follows:
        interpolant = RBFInterpolant(x)
        interpolant.update(f)
        interpolant.evaluate(xi)
    """

    def __init__(self, x, fvals=None,
                 name='thin_plate',
                 func=None, gradient=None,
                 P=lambda x: [1, x[0], x[1]],
                 **params):
        """Data sites (or "centers") with shape (N,dim)
        Default RBF is the thin-plate spline
        Default polynomial is P(x,y) = B0 + B1*x + B2*y
        """
        self.x = x
        self.dist = None
        self.N, self.dim = x.shape
        self.P = P
        self.Q = len(P(x[0,:]))

        def if_not_found(*args,**kwargs):
            print('Warning:',name,'is not a valid RBF name')
        thismodule = sys.modules[self.__module__]
        if func is not None:
            self.func = func
        else:
            self.func = getattr(thismodule,name,if_not_found)
        if gradient is not None:
            self.grad = gradient
        else:
            self.grad = getattr(thismodule,'grad_'+name,if_not_found)

        self.params = params
        self._calculate_LHS()

        if fvals is not None:
            self.update(fvals)
        else:
            self.Ntimes = 0
            self.coefs = None

    def _calculate_LHS(self):
        """Populate the left-hand side of the interpolation system
        for a given set of data sites and specified RBF"""
        Amat = np.zeros((self.N+self.Q, self.N+self.Q))
        for j in range(self.N):
            r2 = np.zeros((self.N,))
            for d in range(self.dim):
                r2 += (self.x[j,d] - self.x[:,d])**2
            Amat[j,:self.N] = self.func(np.sqrt(r2),**self.params)
            pvals = self.P(self.x[j,:])
            Amat[j,self.N:] = pvals
            Amat[self.N:,j] = pvals
        self.A = Amat

    def calculate_separation(self):
        """Calculate radial distances between all combinations of
        centers"""
        self.dist = []
        for j in range(self.N):
            for k in range(j,self.N):
                if j == k:
                    continue
                else:
                    d = self.x[j,:] - self.x[k,:]
                    r = np.sqrt(np.dot(d,d))
                    self.dist.append(r)

    def mean_separation(self):
        if self.dist is None:
            self.calculate_separation()
        return np.mean(self.dist)

    def min_separation(self):
        if self.dist is None:
            self.calculate_separation()
        return np.min(self.dist)

    def max_separation(self):
        if self.dist is None:
            self.calculate_separation()
        return np.max(self.dist)

    def update(self,fvals):
        """Calculates interpolation coefficients for a set of function
        values associated with the data sites. Values should have shape
        (N,Ntimes)"""
        if len(fvals.shape) == 1:
            self.Ntimes = 1
            N = len(fvals)
            fvals = fvals.reshape((N,1))
        else:
            N,self.Ntimes = fvals.shape
        assert(N == self.N)
        coefs = np.zeros((self.N+self.Q,self.Ntimes))
        RHS = np.zeros(np.shape(coefs))
        RHS[:N,:] = fvals
        for itime in range(self.Ntimes):
            coefs[:,itime] = np.linalg.solve(self.A, RHS[:,itime])
        self.coefs = coefs

    def evaluate(self,xi):
        """Evaluates interpolant at xi for all times"""
        assert(self.coefs is not None)
        r2 = np.zeros((self.N,))
        interpvec = np.zeros((self.N+self.Q,))
        for d in range(self.dim):
            r2 += (xi[d] - self.x[:,d])**2
        interpvec[:self.N] = self.func(np.sqrt(r2),**self.params)
        interpvec[self.N:] = self.P(xi)
        return np.dot(interpvec, self.coefs)

    def evaluate_gradient(self,xi,axis):
        """Evaluates gradient of the interpolant at xi for all times, 
        where axisname is 'x', 'y', or 'z'.
        Assumes that the lowest-order polynomial terms are [1, x, y]
        """
        assert(self.coefs is not None)
        assert(self.grad is not None)
        assert(axis < self.dim)
        delta = np.zeros((self.N,self.dim))
        r2 = np.zeros((self.N,))
        interpvec = np.zeros((self.N+self.Q,))
        for d in range(self.dim):
            delta[:,d] = xi[d] - self.x[:,d]
            r2 += delta[:,d]**2
        interpvec[:self.N] = delta[:,axis] * self.grad(np.sqrt(r2),**self.params)
        interpvec[self.N+axis+1] = 1.0
        return np.dot(interpvec, self.coefs)

    def evaluate_field(self,xvec,yvec):
        """Calls evaluate for all points on a regular grid defined by
        x and y vectors"""
        field = np.zeros((len(xvec),len(yvec),self.Ntimes))
        for i,xi in enumerate(xvec):
            for j,yj in enumerate(yvec):
                field[i,j,:] = self.evaluate([xi,yj])
        return field

    def evaluate_gradient_field(self,xvec,yvec):
        """Calls evaluate for all points on a regular grid defined by
        x and y vectors"""
        field = np.zeros((len(xvec),len(yvec),self.Ntimes))
        for i,xi in enumerate(xvec):
            for j,yj in enumerate(yvec):
                field[i,j,:] = self.evaluate_gradient([xi,yj])
        return field


class LOOCV(object):
    """Driver for Leave-One-Out Cross-Validation (LOOCV), which uses
    the RBFInterpolant class.
    """

    def __init__(self,x,**kwargs):
        self.x = x
        self.N, dim = x.shape
        self.xsubsets = [ np.concatenate((x[:i],x[i+1:]))
                          for i in range(self.N) ]
        self.interpolants = [ RBFInterpolant(xsub, **kwargs)
                              for xsub in self.xsubsets ]
        self.fsubsets = None

    def calculate(self,fvals):
        N,self.Ntimes = fvals.shape
        assert(N == self.N)
        self.fref = fvals
        self.fsubsets = [ np.concatenate((fvals[:i,:],fvals[i+1:,:]),
                                         axis=0)
                          for i in range(self.N) ]
        for interp,f in zip(self.interpolants,self.fsubsets):
            interp.update(f)
        self.finterp = [ 
            interp.evaluate(xi)
            for xi,interp in zip(self.x,self.interpolants)
        ]
        self.errors = np.array(self.finterp) - self.fref
        return self.errors
