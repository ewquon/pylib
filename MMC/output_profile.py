# This produces standardized ASCII output for NCAR data assessment in
# the DOE Mesoscale to Microscale Modeling project.
#
# Original FORTRAN code by Branko Kosovic (NCAR) 06-03-2015
# Python code by Eliot Quon (NREL) 08-01-2018

from __future__ import print_function
import os
import io # for buffering
#print('Default buffer size: ',io.DEFAULT_BUFFER_SIZE)

import numpy as np
import pandas as pd



profile_data = [
    'u',      # streamwise velocity compnent, units [m/s]
    'v',      # cross-stream velocity compnent, units [m/s]
    'w',      # vertical velocity compnent, units [m/s]
    'th',     # potential or virtual potential temperature, units [K]
    'p',      # pressure [mbar]
    'tke',    # turbulent kinetic energy, units [m^2/s^2] 
    'tau11',  # Reynolds or subgrid/subfilter stress component 11 <uu>, units [m^2/s^2] 
    'tau12',  # Reynolds or subgrid/subfilter stress component 12 <uv>, units [m^2/s^2] 
    'tau13',  # Reynolds or subgrid/subfilter stress component 13 <uw>, units [m^2/s^2] 
    'tau22',  # Reynolds or subgrid/subfilter stress component 22 <vv>, units [m^2/s^2] 
    'tau23',  # Reynolds or subgrid/subfilter stress component 23 <vw>, units [m^2/s^2] 
    'tau33',  # Reynolds or subgrid/subfilter stress component 33 <ww>, units [m^2/s^2] 
    'hflux',  # Reynolds or subgrid/subfilter vertical heat flux <wth>, units [Km/s] 
]

surface_data = [
    'ustar',  # surface friction velocity, units [m/s]
    'z0',     # surface roughness, units [m]
    'T0',     # skin temperature, units [K]
    'qwall',  # surface flux, units [Km/s]
]


header = """INSTITUTION:{institution:s}
   LOCATION:{location:s}
   LATITUDE:{lat:10.4f}
  LONGITUDE:{long:10.4f}
   CODENAME:{codename:s}
   CODETYPE:{codetype:s}
   CASENAME:{casename:s}
  BENCHMARK:{benchmark:s}
     LEVELS:{levels:7d}
"""

snapshot = """
       DATE:{date:s}
       TIME:  {time:s}
FRICTION VELOCITY [m/s] = {ustar:10.5f}
SURFACE ROUGHNESS [m]   = {z0:10.5f}
SKIN TEMPERATURE [K]    = {T0:10.5f}
SURFACE FLUX [Km/s]     = {qwall:10.5f}
             Z [m]           U [m/s]           V [m/s]           W [m/s]            TH [K]           P [mbar]    TKE [m^2/s^2]   TAU11 [m^2/s^2]   TAU12 [m^2/s^2]   TAU13 [m^2/s^2]   TAU22 [m^2/s^2]   TAU23 [m^2/s^2]   TAU33 [m^2/s^2]      HFLUX [Km/s]
"""
row = 4*'{:18.3f}' + 2*'{:18.2f}' + '{:18.3f}' + 7*'{:18.5f}' + '\n'


class MMCdata(object):

    def __init__(self, profile, surface, info, na_values=-999.0):
        """Create object for standardized MMC data output

        profile: pandas.DataFrame
            Dataframe with datetime index and columns corresponding to
            the data arrays described above; a 'z' column describes the
            height AGL. 
        surface: pandas.DataFrame
            Dataframe with datetime index and columns corresponding to
            the surface data arrays described above.
        info: dict
            Dictionary with header information, including institution,
            location, codename, casetype, casename, latitude, and 
            longitude.
        na_values: float
            Default value if data do not exist.
        """
        if surface is not None:
            self.t = surface.index
            assert(np.all(self.t == profile.index.unique()))
            self.surface = surface
        else:
            self.t = profile.index.unique()
            self.surface = pd.DataFrame(index=self.t)
        self.z = profile['z'].unique()
        self.info = info
        self.N = len(self.z)
        self.info['levels'] = self.N
        self.na = na_values

        self.profile = profile
        self.profile.index.name = 'date_time'
        self.profile = self.profile.reset_index()
        self.profile = self.profile.sort_values(['date_time','z'])
        self.profile = self.profile.set_index('date_time')

        # fill in na_values for unavailable data
        for name in surface_data:
            if name not in self.surface.columns:
                self.surface[name] = self.na
        self.surface = self.surface.fillna(self.na)

        for name in profile_data:
            if name not in self.profile.columns:
                self.profile[name] = self.na
        self.profile = self.profile.fillna(self.na)


    def write_ascii(self,fname,surface=None,profile=None):
        """Write out ascii data (for all available dates by default)"""
        if surface is None:
            surface = self.surface
        if profile is None:
            profile = self.profile
        # convert strings to uppercase
        infodata = self.info.copy()
        for key,val in infodata.items():
            if isinstance(val,str):
                infodata[key] = val.upper()
        #with open(fname,'w') as f:
        with io.open(fname,'w') as f:
            f.write(header.format(**infodata))
# df.to_dict(orient='records') returns a list of dicts, with one dict
# per row of the dataframe
#-- THIS IS SLOW
#            for t_i, hdr in zip(self.t,
#                                     self.surface.to_dict(orient='records')):
#-- THIS IS STILL SLOW
#            for t_i, rowdata in self.surface.iterrows():
            for i, (t_i, rowdata) in enumerate(surface.iterrows()):
                hdr = rowdata.to_dict()
                hdr['date'] = str(t_i.date())  # format [YYYY-MM-DD]
                hdr['time'] = str(t_i.time())  # format [HH:MM:SS] UTC
                f.write(snapshot.format(**hdr))
#-- MAJOR SLOWDOWN HERE
#                df = self.profile.loc[self.profile.index == t_i]
                df = profile.iloc[i*self.N:(i+1)*self.N]
                for z_i, zdata in df.iterrows():
                    zdata = zdata[['z']+profile_data].values
                    f.write(row.format(*zdata))
        print('Wrote '+fname)


    def write_ascii_days(self,dpath,fname_fmt='%Y%m%d.dat'):
        """Call write_ascii for individual days"""
        days = np.unique(self.t.date)
        print('Days in data:',[str(day) for day in days])
        surface_date = self.surface.index.date
        profile_date = self.profile.index.date
        for day in days:
            surf = self.surface.loc[surface_date == day]
            prof = self.profile.loc[profile_date == day]
            fname = day.strftime(fname_fmt)
            #print('Writing out '+fname)
            fpath = os.path.join(dpath,fname)
            if os.path.isfile(fpath):
                print('Skipping existing file',fname)
                continue
            self.write_ascii(fpath, surface=surf, profile=prof)

