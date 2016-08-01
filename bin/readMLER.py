#!/usr/local/bin/python
#
# Utility for processing MLER-generated coefficients files
#
# Eliot Quon (eliot.quon@nrel.gov)
# Written 6/9/16
#
import sys, os
import numpy as np
import matplotlib.pyplot as plt

import argparse # {{{
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, \
        description="""Reads and plots MLER-generated focused wave
Expected superposition form:

  waveElev = sum( sqrt(2*SpectAmp*dw) * cos( -k*(x-X0) + w*(t-T0) + Phase) )

Note that phases are specified relative to X0 and T0 in space and time.""")
parser.add_argument('fname', metavar='filename', 
        type=str, default='',
        help='path to MLER coefficients file')
parser.add_argument('-v', '--verbose', action='store_const',
        const=True, default=True)

outputArgs = parser.add_argument_group('output arguments')
outputArgs.add_argument('-a', '--all', action='store_const',
        const=True, default=False,
        help='make all plots')
outputArgs.add_argument('-i', '--inletTrace', action='store_const',
        const=True, default=False,
        help='plot z vs t at inlet')
outputArgs.add_argument('-d', '--deviceTrace', action='store_const',
        const=True, default=False,
        help='plot z vs t at device location')
outputArgs.add_argument('-0', '--initialProfile', action='store_const',
        const=True, default=False,
        help='plot z vs x at t=0')
outputArgs.add_argument('-p', '--peakProfile', action='store_const',
        const=True, default=False,
        help='plot z vs x at t=tpeak-T0')
outputArgs.add_argument('--movie', action='store_const',
        const=True, default=False,
        help='generate image snapshots for assembly into a movie')
#outputArgs.add_argument('--min', action='store_const',
#        const=True, default=False,
#        help='calculate the minimum wave elevation')

constArgs = parser.add_argument_group('constants (may be overwritten by input file header)')
constArgs.add_argument('-t','--tpeak', metavar='tp', 
        type=float, default=150.0,
        help='time (simulated) of peak response [s]  (default=150 s)')
constArgs.add_argument('--t0', metavar='T0', 
        type=float, default=0.0,
        help='time (MLER) of peak response [s]  (default=0 s)')
constArgs.add_argument('--x0', metavar='X0', 
        type=float, default=0.0,
        help='peak (device) position [m]  (default=0 m)')
constArgs.add_argument('--xin', metavar='xi',
        type=float, default=-600.0,
        help='domain inlet location [m]  (default=-600 m)')
constArgs.add_argument('--xout', metavar='xo',
        type=float, default=600.0,
        help='domain outlet location [m] corresponding to the mean water line reference location  (default= 600 m)')

args = vars(parser.parse_args())
fname = args['fname']
if not os.path.isfile(fname): sys.exit('MLER coefficients file {} not found'.format(fname))

plot_z_vs_t_at_inlet  = args['inletTrace']
plot_z_vs_t_at_device = args['deviceTrace']
plot_z_vs_x_at_tinit  = args['initialProfile']
plot_z_vs_x_at_tpeak  = args['peakProfile']
if args['all']:
    plot_z_vs_t_at_inlet = plot_z_vs_t_at_device = plot_z_vs_x_at_tinit = plot_z_vs_x_at_tpeak = True
generate_movie_snapshots = args['movie']
#calc_eta_min = args['min']

x0 = args['x0']
t0 = args['t0']
tpeak = args['tpeak']
xin = args['xin']
xout = args['xout']

verbose = args['verbose']
# }}}
#--------------------------------------------
# hard-coded inputs
Nt = 501 # plot resolution
Nx = 501
moviedir = 'snapshots'
#--------------------------------------------

trange = np.abs(tpeak)

# read w, A, p, k from MLER input file# {{{
w = [] # frequency, rad/s
S = [] # spectral amplitude, m^2
p = [] # phase, rad
k = [] # wavenumber, rad^2/m
with open(fname,'r') as f:
    for line in f:
        if line.strip()=='': break
        if line.startswith('#'):
            if verbose: sys.stdout.write(line)
            line = line.split()
            if not line[0]=='#': continue
            try:
                param = line[1]
                if param.startswith('dW'):
                    dw = float(line[2])
                    print '<<<<<< Setting',param,dw,'>>>>>>'
                elif param.startswith('T0'):
                    newval = float(line[2])
                    if not t0==newval:
                        t0 = newval
                        print '<<<<<< Overriding',param,t0,'>>>>>>'
                elif param.startswith('X0'):
                    newval = float(line[2])
                    if not x0==newval:
                        x0 = newval
                        print '<<<<<< Overriding',param,x0,'>>>>>>'
            except: pass
            continue

        wi,Si,pp,ki = [ float(val) for val in line.split() ]
        if Si==0: continue
        w.append(wi)
        S.append(Si)
        p.append(pp)
        k.append(ki)

N = len(w)
w = np.array(w)
if not 'dw' in vars(): dw = w[1] - w[0]
print 'dw =',dw
S = np.array(S)
A = np.sqrt(2*S*dw)  #np.array(S)**0.5
p = np.array(p)
k = np.array(k)
T = 2*np.pi / w
L = 2*np.pi / k
c = w/k# }}}

print N,'subwaves read'
print 'omega:      [',np.min(w),',',np.max(w),']'
print 'amplitude:  [',np.min(A),',',np.max(A),']'
print 'phase:      [',np.min(p),',',np.max(p),']'
print 'period:     [',np.min(T),',',np.max(T),']'
print 'wavelength: [',np.min(L),',',np.max(L),']'
print 'wavespeed:  [',np.min(c),',',np.max(c),']  at w=',w[np.argmax(c)],'rad/s'

imax = np.argmax(A)
print '\nspectral peak at w={:f} rad/s :  T={:f} s, A={:f} m, L={:f} m, phi={:f} rad\n'.format(w[imax],T[imax],A[imax],L[imax],p[imax])

eta_max = np.sum( A*np.cos(p) )
print '(max) wave elevation at t=0 :',eta_max

dzdx_max = np.sum( A*k*np.sin(p) )
print '(max) wave slope at t=0 :',dzdx_max

#
# Final MLER superposition form (06/2016):
#   waveElev = sum( sqrt(2*SpectAmp*dw) * cos( -k*(x-X0) + w*(t_MLER-T0) + Phase) )
#
# Star superposition form:
#   z_star = sum( A * cos( k*(x-xref) - w*t_star - phi_star + pi/2 ) )
# where in general t_star = t_MLER + t_peak
#
#def inputsubwave(x,t,phi): return A*np.cos( k*x - w*t - phi ) # new formulation with different sign for phi (2/10/16)
def inputsubwave(x,t,phi): return A*np.cos( -k*(x-x0) + w*(t-t0) + phi ) # (updated 6/8/16)
def Star_subwave(x,t,phi): return A*np.sin( k*(x-xout) - w*t + phi + np.pi ) # for ref pt at (xout,0,0) (2/24/16)

# DEBUG
#x = np.linspace(xin,xout,Nx)
#plt.plot( x, A[imax]*np.cos( -k[imax]*(x-x0) + w[imax]*(-tpeak) + p[imax] ) ) # input subwave at t=t0-tpeak
#plt.grid()

# phase correction for Star
#centershift = -k*xout + np.pi/2
#center_timeshift = centershift - w*tpeak
center_timeshift =  k*xout + w*tpeak - np.pi/2 # updated 6/8/16
center_timeshift2 = k*(xout-x0) + w*(t0+tpeak) - 2*p - 3*np.pi/2 # updated 6/8/16, NEED TO FLIP SIGN ON A
center_timeshift3 = k*(xout-x0) + w*(t0+tpeak) - np.pi/2 # 6/8/16, no sign-flip on A, need phi=-p+center_timeshift3

# TODO: This doesn't necessarily give the correct global minima
#if calc_eta_min:# {{{
#    from scipy.optimize import minimize
#
#    deriv_at_x0 = lambda t: np.sum( -A*w*np.sin( w*(t-t0) + p ) )
#    tguess = np.linspace(-trange,trange,10*Nt)
#    devicetrace = [ np.sum( inputsubwave(x0,ti,p) ) for ti in tguess ]
#    t_min_guess = tguess[np.argmin(devicetrace)]
#    result = minimize( deriv_at_x0, t_min_guess )
#    try:
#        print result
#        t_min = result.x
#        eta_min_x0 = np.sum( inputsubwave(x0,t_min,p) )
#        print 'min wave elevation at x=0 :',eta_min_x0,'at t=',t_min,'s'
#    except:
#        print result
#        print 'scalar minimization at x=0 failed'
#
#    deriv_at_t0 = lambda x: np.sum( A*k*np.sin( -k*(x-x0) + p ) )
#    xguess = np.linspace(xin,xout,10*Nx)
#    peakprofile = [ np.sum( inputsubwave(xi,t0,p) ) for xi in xguess ]
#    x_min_guess = xguess[np.argmin(peakprofile)]
#    result = minimize( deriv_at_t0, x_min_guess )
#    try:
#        print result
#        x_min = result.x
#        eta_min_t0 = np.sum( inputsubwave(x_min,t0,p) )
#        print 'min wave elevation at t=0 :',eta_min_t0,'at x=',x_min,'m'
#    except:
#        print result
#        print 'scalar minimization at t=0 failed'# }}}

#
# plot
#

if plot_z_vs_t_at_inlet:# {{{

    t_in        = np.linspace(-trange,trange,Nt)
    t           = np.linspace(0,2*trange,Nt)
    zinput0     = np.zeros((Nt))
    zinput1     = np.zeros((Nt))
    z_star0     = np.zeros((Nt))
    z_star1     = np.zeros((Nt))
    z_star2     = np.zeros((Nt))
    for i in range(Nt):
        zinput0[i] =  np.sum( inputsubwave( xin, t_in[i],     p ) )
        zinput1[i] =  np.sum( inputsubwave( xin, t[i]-tpeak,  p ) )
        z_star0[i] =  np.sum( Star_subwave( xin, t[i],        p+center_timeshift ) )
        z_star1[i] = -np.sum( Star_subwave( xin, t[i],        p+center_timeshift2 ) )
        z_star2[i] =  np.sum( Star_subwave( xin, t[i],       -p+center_timeshift3 ) )
    fig, (ax0,ax1) = plt.subplots(2)
    ax0.plot(t_in,zinput0,'k-',label='input')
    ax0.set_ylabel('INPUT z(x={:.1f},t)'.format(xin))
    ax0.legend(loc='best')
    ax1.plot(t,zinput1,'k',linewidth=2,label='input (timeshift)')
    ax1.plot(t,z_star0,label='Star (center+timeshift)')
    ax1.plot(t,z_star1,label='Star (center+timeshift, -A)')
    ax1.plot(t,z_star2,label='Star (center+timeshift, A, -p)')
    ax1.set_xlabel('TIME t')
    ax1.set_ylabel('SIMULATED HEIGHT z(x={:.1f},t)'.format(xin))
    ax1.legend(loc='best')
    fig.suptitle('Profile height at INLET')# }}}

if plot_z_vs_t_at_device:# {{{

    t_in        = np.linspace(-trange,trange,Nt)
    t           = np.linspace(0,2*trange,Nt)
    zinput0     = np.zeros((Nt))
    zinput1     = np.zeros((Nt))
    z_star0     = np.zeros((Nt))
    z_star1     = np.zeros((Nt))
    z_star2     = np.zeros((Nt))
    for i in range(Nt):
        zinput0[i] =  np.sum( inputsubwave( x0, t_in[i],     p ) )
        zinput1[i] =  np.sum( inputsubwave( x0, t[i]-tpeak,  p ) )
        z_star0[i] =  np.sum( Star_subwave( x0, t[i],        p+center_timeshift ) )
        z_star1[i] = -np.sum( Star_subwave( x0, t[i],        p+center_timeshift2 ) )
        z_star2[i] =  np.sum( Star_subwave( x0, t[i],       -p+center_timeshift3 ) )
    fig, (ax0,ax1) = plt.subplots(2)
    ax0.plot(t_in,zinput0,'k-',label='input')
    ax0.set_ylabel('INPUT z(x={:.1f},t)'.format(x0))
    ax0.legend(loc='best')
    ax1.plot(t,zinput1,'k',linewidth=2,label='input (timeshift)')
    ax1.plot(t,z_star0,label='Star (center+timeshift)')
    ax1.plot(t,z_star1,label='Star (center+timeshift, -A)')
    ax1.plot(t,z_star2,label='Star (center+timeshift, A, -p)')
    ax1.set_xlabel('TIME t')
    ax1.set_ylabel('SIMULATED HEIGHT z(x={:.1f},t)'.format(x0))
    ax1.legend(loc='best')
    fig.suptitle('Profile height at DEVICE')# }}}

if plot_z_vs_x_at_tinit:# {{{

    x = np.linspace(xin,xout,Nx)
    zinput0 = np.zeros((Nx))
    zinput1 = np.zeros((Nx))
    z_star0 = np.zeros((Nx))
    z_star1 = np.zeros((Nx))
    z_star2 = np.zeros((Nx))
    for i in range(Nx):
        zinput0[i] =  np.sum( inputsubwave( x[i], t0,        p ) )
        zinput1[i] =  np.sum( inputsubwave( x[i], t0-tpeak,  p ) )
        z_star0[i] =  np.sum( Star_subwave( x[i], 0.,        p+center_timeshift ) )
        z_star1[i] = -np.sum( Star_subwave( x[i], 0.,        p+center_timeshift2 ) )
        z_star2[i] =  np.sum( Star_subwave( x[i], 0.,       -p+center_timeshift3 ) )
    fig, (ax0,ax1) = plt.subplots(2,sharex=True)
    ax0.plot(x,zinput0,'b-',label='input (t0={:.1f})'.format(t0))
    ax0.plot(x,zinput1,'k-',linewidth=2,label='input (timeshift)')
    ax0.set_ylabel('INPUT z(x,t=0)')
    ax0.legend(loc='best')
    ax1.plot(x,zinput1,'k-',linewidth=2,label='input (timeshift)')
    ax1.plot(x,z_star0,label='Star (center+timeshift)')
    ax1.plot(x,z_star1,label='Star (center+timeshift, -A)')
    ax1.plot(x,z_star2,label='Star (center+timeshift, A, -p)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('SIMULATED HEIGHT z(x,t=0)')
    titlestr = 'INITIAL profile height'
    ax1.legend(loc='best')
    fig.suptitle(titlestr)# }}}

if plot_z_vs_x_at_tpeak:# {{{

    x = np.linspace(xin,xout,Nx)
    zinput0 = np.zeros((Nx))
    z_star0 = np.zeros((Nx))
    z_star1 = np.zeros((Nx))
    z_star2 = np.zeros((Nx))
    for i in range(Nx):
        zinput0[i] =  np.sum( inputsubwave( x[i], t0,     p ) )
        z_star0[i] =  np.sum( Star_subwave( x[i], tpeak,  p+center_timeshift ) )
        z_star1[i] = -np.sum( Star_subwave( x[i], tpeak,  p+center_timeshift2 ) )
        z_star2[i] =  np.sum( Star_subwave( x[i], tpeak, -p+center_timeshift3 ) )
    fig, (ax0,ax1) = plt.subplots(2,sharex=True)
    ax0.plot(x,zinput0,'k-',label='input')
    ax0.set_ylabel('INPUT z(x,t={:.1f})'.format(tpeak))
    ax0.legend(loc='best')
    ax1.plot(x,zinput0,'k-',linewidth=2,label='input')
    ax1.plot(x,z_star0,label='Star (center+timeshift)')
    ax1.plot(x,z_star1,label='Star (center+timeshift, -A)')
    ax1.plot(x,z_star2,label='Star (center+timeshift, A, -p)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('SIMULATED HEIGHT z(x,t={:.1f})'.format(tpeak))
    ax1.legend(loc='best')
    fig.suptitle('PEAK profile height')# }}}

#---------------------
plt.show()

# do this last...
if generate_movie_snapshots:# {{{
    x = np.linspace(xin,xout,Nx)
    if generate_movie_snapshots:
        ymax = max(np.ceil(eta_max),11)
        plt.figure()
        plt.hold(False)
        name = '.'.join(fname.split('.')[:-1])
        if not os.path.isdir(moviedir):
            print 'Creating output dir :',moviedir
            os.makedirs(moviedir)
        if verbose: print '\nSaving movie snapshots to',moviedir
    for t in np.arange(-trange,trange+1):
        eta = np.zeros(x.shape)
        for i in range(len(x)):
            eta[i] = np.sum( inputsubwave(x[i],t,p) )
        plt.plot(x, eta)
        plt.ylim((-ymax,ymax))
        plt.title('t=%f'%t)
        imgname = moviedir + os.sep + '{:s}_{:04d}.png'.format(name,int(t+trange))
        plt.savefig(imgname)
        print '  wrote',imgname
    print 'To generate an mp4:'
    print '  ffmpeg -f image2 -r 5 -i snapshots/coeffs_%04d.png -vcodec mpeg4 -y movie.mp4'
# }}}

