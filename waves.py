#!/usr/bin/python
import sys
import numpy as np

# TODO: integrate fenton1985 with this module
sys.path.append('/Users/equon/WEC/fifthOrderWave')
import fenton1985

g = 9.81
theories = ['linear','fenton1985']

# Numerical parameters
DEBUG = False
MAXITER = 25
TOL = 1e-14

def solve_k(w,d=70.0,guess=None,H=-1,theory='linear'):
    """ Solve for the wavenumber (k) satisfying the linear dispersion relation at a given depth (d) and frequency (w).
    An optional guess may be provided, otherwise the initial k value will be estimated from the deep-water limit.
    """
    if not theory in theories: 
        print 'Available theories:',theories
        return -1

    if guess: ki = guess
    else: ki = w*w/g

    niter = 0
    if theory=='linear':
        if DEBUG: print 'Solving k for w={:f} rad/s at d={:f} m'.format(w,d)
        F,J = _linear_dispersion_newton_step(w,ki,d)
        relax = 1.0
    elif theory=='fenton1985':
        if DEBUG: print 'Solving k for w={:f} rad/s, H={:f} m at d={:f} m'.format(w,H,d)
        print 'WARNING this is deprecated, use scipy.optimize instead'
        F,J = _fenton1985_dispersion_step(w,ki,H,d)
        relax = 0.5
    if DEBUG: print '  iter',niter,'resid',F,' k=',ki

    while niter < MAXITER and np.abs(F) > TOL:

        #ki += -F/J
        ki += -F/J * relax
        niter += 1
        if theory=='linear':
            F,J = _linear_dispersion_newton_step(w,ki,d)
        elif theory=='fenton1985':
            F,J = _fenton1985_dispersion_step(w,ki,H,d)

        if DEBUG: print '  iter',niter,'resid',F,' k=',ki

    if niter >= MAXITER: print 'WARNING: max newton iterations reached!'
    if DEBUG: 
        if niter < MAXITER: print '  Newton iteration converged in',niter,'steps'

    return ki

def solve_fifthorder_k(w,H,d=70.0,guess=None):
    #solve_k(w,d,guess,H=H,theory='fenton1985')
    from scipy.optimize import fsolve
    def eqn23(k): # TODO: handle mean current speed not 0
        C0,C2,C4 = fenton1985.evalC(k*d)
        return -w/(g*k)**0.5 + C0 + (k*H/2)**2*C2 + (k*H/2)**4*C4
    if not guess: guess = w**2/g # deep water approximation, use as a starting guess
    k = fsolve(eqn23,guess)
    if isinstance(k,np.ndarray): k = k[0] # depending on version, may return array or scalar

    return k


#---------------------------------------------------------------------------------------------------
# 
# iterative solution step for various theories when applying Newton's method
#

def _linear_dispersion_newton_step(w,k,d):
    kd = k*d
    tkd = np.tanh(kd)
    return k*tkd - w**2/g, kd*(1-tkd**2) + tkd # F(k), dF/dk

def _fenton1985_dispersion_step(w,k,H,d):
    """ Solve eqn 23 in Fenton1985, assuming no mean current (i.e., c_E = 0) using Newton's method.
    This approach was found to be numerically unstable unless a relaxation factor was applied at each iteration; however this required > 100 iterations for convergence in the cases tested. 
    """
    kd = k*d
    S = 1./np.cosh(2*kd)
    C0 = np.tanh(kd)**0.5
    dC0 = d/(2*C0*np.cosh(kd)**2)
    C2expr = (2.+7.*S**2)/(4.*(1.-S)**2)
    C2 = C0 * C2expr
    t2kd = np.tanh(2*kd)
    dC2 = \
            dC0 * C2expr \
            + C0 *( -7.*d*t2kd*S**2/(1.-S)**2 - d*t2kd*(7*S**2+2.)*S/(1.-S)**3 )
    C4expr = (4. + 32.*S - 116.*S**2 - 400.*S**3 - 71.*S**4 + 146.*S**5) / (32*(1.-S)**5)
    C4 = C0 * C4expr
    dC4 = \
            dC0 * C4expr \
            + C0 * ( -1460*d*t2kd*S**5 + 568*d*t2kd*S**4 \
                + 2400*d*t2kd*S**3 + 464*d*t2kd*S**2 - 64*d*t2kd*S ) / (32*(1.-S)**5) \
            - C0 * ( 5*d*t2kd*S*(146*S**5-71*S**4-400*S**3-116*S**2+32*S+4) ) \
                    / (16*(1.-S)**6)
    kH_2 = k*H/2.
    F = -w/(g*k)**0.5 + C0 + kH_2**2*C2 + kH_2**4*C4
    J = dC0 + H*kH_2*C2 + kH_2**2*dC2 + 2*H*kH_2**3*C4 + kH_2**4*dC4

    if DEBUG:
        print ' C0:',C0
        print ' C2:',C2
        print ' C4:',C4
        print 'dC0:',dC0
        print 'dC2:',dC2
        print 'dC4:',dC4
    return F,J

