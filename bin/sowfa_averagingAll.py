#!/usr/bin/env python
import sys
from SOWFA.postProcessing.averaging import read
import refs
import matplotlib.pyplot as plt
import numpy as np

heights = [90.]
tavg_window = 600.
dt = 1.0
SFS = True

ref = refs.read()

data = read( *sys.argv[1:] )

#------------------------------------------------------------------------------
tavg,TIx_hist,TIy_hist,TIz_hist,TIdir_hist,TIxyz_hist,TKE_hist = data.calcTI_hist( heights, tavg_window, dt, SFS )
for ih,z in enumerate(heights):
    plt.plot( tavg, TIdir_hist[:,ih], label='z={:.1f}'.format(z) )
    fname = 'TIhist_z{:.1f}.csv'.format(z)
    np.savetxt( fname, np.vstack((tavg,TIdir_hist[:,ih])).T, delimiter=',', header='Time,TI' )
    print 'wrote',fname

plt.legend(loc='best')
plt.savefig('TIhist.png')

#==========
plt.show()
#==========

data.calcTI( heights )
for i,h in enumerate(heights):
    print 'z=',h,':'
    print '  TIx/y/z =',data.TIx[i],data.TIy[i],data.TIz[i]
    print '  TIdir   =',data.TIdir[i]
    print '  TIxyz   =',data.TIxyz[i],' ( TKE=',data.TKE[i],')'


#------------------------------------------------------------------------------
fig00,ax00 = plt.subplots(ncols=2)

ax00[0].plot( data.U_mean[-1,:], data.hLevelsCell, label=r'$U$' )
ax00[0].plot( data.V_mean[-1,:], data.hLevelsCell, label=r'$V$' )
ax00[0].plot( data.W_mean[-1,:], data.hLevelsCell, label=r'$W$' )
if ref:
    try:
        z = data.hLevelsCell
        refprofile = ref.Uref * (z/ref.zref)**ref.shear
        ax00[0].plot( refprofile, z, 'k--', lw=2, label='input' )
    except (NameError,TypeError): pass
ax00[0].set_xlabel('Velocity (m/s)')
ax00[0].set_ylabel('Height (m)')
ax00[0].legend(loc='best')

ax00[1].plot( data.T_mean[-1,:], data.hLevelsCell )
ax00[1].set_xlabel('Temperature (K)')

fig00.suptitle('Resolved Mean Quantities')
fig00.savefig('Profiles_Mean.png')

#------------------------------------------------------------------------------
fig0,ax0 = plt.subplots(ncols=2)

windMag = np.sqrt( data.U_mean[-1,:]**2 + data.V_mean[-1,:]**2 )
windDir = np.arctan2( -data.U_mean[-1,:], -data.V_mean[-1,:] )*180.0/np.pi

ax0[0].plot( windMag, data.hLevelsCell, label=r'$U$' )
if ref:
    try:
        z = data.hLevelsCell
        refprofile = ref.Uref * (z/ref.zref)**ref.shear
        ax0[0].plot( refprofile, z, 'k--', lw=2, label='input' )
    except (NameError,TypeError): pass
ax0[0].set_xlabel('Horizontal Velocity (m/s)')
ax0[0].set_ylabel('Height (m)')
#ax0[0].legend(loc='best')

ax0[1].plot( windDir, data.hLevelsCell )
ax0[1].set_xlabel('Direction (deg)')

if ref:
    try:
        z0 = ref.zref - ref.D/2
        z1 = ref.zref + ref.D/2
        dir0 = np.interp(z0,data.hLevelsCell,windDir)
        dir1 = np.interp(z1,data.hLevelsCell,windDir)
        dirHH = np.interp(ref.zref,data.hLevelsCell,windDir)
        ax0[1].axhline(z0,linestyle='-',color='0.5')
        ax0[1].axhline(z1,linestyle='-',color='0.5')
        ax0[1].annotate(r' $\Delta = {:.1f}^\circ$'.format(np.abs(dir1-dir0)),(dirHH,ref.zref))
    except NameError: pass

fig0.suptitle('Resolved Mean Quantities')
fig0.savefig('Profiles_MeanUdir.png')

#------------------------------------------------------------------------------
fig1,ax1 = plt.subplots(ncols=2)

ax1[0].plot( data.uu_mean[-1,:], data.hLevelsCell, label=r"$<u'u'>$" )
ax1[0].plot( data.vv_mean[-1,:], data.hLevelsCell, label=r"$<v'v'>$" )
ax1[0].plot( data.ww_mean[-1,:], data.hLevelsCell, label=r"$<w'w'>$" )
ax1[0].set_xlabel('Normal stresses (m^2/s^2)')
ax1[0].set_ylabel('Height (m)')
ax1[0].legend(loc='best')

ax1[1].plot( data.uv_mean[-1,:], data.hLevelsCell, label=r"$<u'v'>$" )
ax1[1].plot( data.uw_mean[-1,:], data.hLevelsCell, label=r"$<u'w'>$" )
ax1[1].plot( data.vw_mean[-1,:], data.hLevelsCell, label=r"$<v'w'>$" )
ax1[1].set_xlabel('Shear stresses (m^2/s^2)')
ax1[1].legend(loc='best')

fig1.suptitle('Resolved Fluctuating Quantities')
fig1.savefig('Profiles_Fluc.png')

#------------------------------------------------------------------------------
fig2,ax2 = plt.subplots(ncols=2)

ax2[0].plot( data.R11_mean[-1,:], data.hLevelsCell, label=r"$R_{11}$" )
ax2[0].plot( data.R22_mean[-1,:], data.hLevelsCell, label=r"$R_{22}$" )
ax2[0].plot( data.R33_mean[-1,:], data.hLevelsCell, label=r"$R_{33}$" )
ax2[0].set_xlabel('Normal stresses (m^2/s^2)')
ax2[0].set_ylabel('Height (m)')
ax2[0].legend(loc='best')

ax2[1].plot( data.R12_mean[-1,:], data.hLevelsCell, label=r"$R_{12}$" )
ax2[1].plot( data.R13_mean[-1,:], data.hLevelsCell, label=r"$R_{13}$" )
ax2[1].plot( data.R23_mean[-1,:], data.hLevelsCell, label=r"$R_{23}$" )
ax2[1].set_xlabel('Shear stresses (m^2/s^2)')
ax2[1].legend(loc='best')

fig2.suptitle('Sub-Filter Scale Quantities')
fig2.savefig('Profiles_SFS.png')

#==========
plt.show()
#==========

fname = 'averagingProfiles.csv'
np.savetxt( fname,
        np.vstack((
            data.hLevelsCell,
            data.U_mean[-1,:],
            data.V_mean[-1,:],
            data.W_mean[-1,:],
            data.T_mean[-1,:],
            data.uu_mean[-1,:],
            data.vv_mean[-1,:],
            data.ww_mean[-1,:],
            data.uv_mean[-1,:],
            data.uw_mean[-1,:],
            data.vw_mean[-1,:],
            data.R11_mean[-1,:],
            data.R22_mean[-1,:],
            data.R33_mean[-1,:],
            data.R12_mean[-1,:],
            data.R13_mean[-1,:],
            data.R23_mean[-1,:]
            )).T,
        header='z,U,V,W,T,uu,vv,ww,uv,uw,vw,R11,R22,R33,R12,R13,R23',
        delimiter=',' )
print 'wrote',fname
