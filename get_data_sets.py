# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 23:19:06 2017
@author: R. A. Real-Rangel (Institute of Engineering UNAM, Mexico)
"""

# http://unidata.github.io/netcdf4-python/
from netCDF4 import Dataset, num2date
import datetime as dt
import numpy as np
import pandas as pd


def get_cs_ts(input_file):
    """
    Extracts the data from a direct observation (dobs) files.
    """
    # Import CSV file.
    with open(input_file, 'r') as f:
        textin = f.read()
    textin = textin.replace('# ', '')
#    header = textin[:findchar(textin, '\n')[3]]
    header = textin[:[i for i, ltr in enumerate(textin) if ltr == '\n'][3]]
    varname = header.splitlines()[2].split(',')
    dtype = np.empty(shape=np.shape(varname), dtype='|S6')

    for n, name in enumerate(varname):
        if (name == 'YEAR') or (name == 'MONTH') or (name == 'DAY'):
            dtype[n] = 'int'
        else:
            dtype[n] = 'f4'

    data = np.genfromtxt(
            fname=input_file,
            delimiter=',',
            skip_header=4,
            dtype=dtype,
            names=varname)

    # Extract time series (ts).
    dset_t = np.asarray(
            [dt.datetime(
                    year=i['YEAR'],
                    month=i['MONTH'],
                    day=i['DAY'],
                    hour=8) for i in data])
    dset_val = np.array([i['PRECIP'] for i in data])
#    dset_val[dset_val < 0.1] = 0

    # Generate pandas dataframe
    dobs_ts = pd.DataFrame()
    dobs_ts[input_file.split('/')[-1].split('.')[0]] = dset_val
    dobs_ts.index = dset_t
    dobs_ts = dobs_ts[~dobs_ts.index.duplicated(keep='first')]
    return(dobs_ts)


#def get_da_ts(da_db_dir):
#    da_file_list = [os.path.join(da_db_dir, f)
#                    for f in os.listdir(da_db_dir)
#                    if f.endswith('.nc')]
#    
#    for daf, da_file in enumerate(da_file_list):
#        with Dataset(da_file, 'r') as f:
#            lat = f['latitude'][:]
#            lon = f['longitude'][:]
#
#            if 'ts_x' in locals() or 'ts_x' in globals():
#                ts_x = np.vstack((ts_x, f['PRECTOTLAND'][:]))
#                ts_t = np.append(ts_t, num2date(f['time'][:],
#                                      units='hours since 1980-01-01 00:00',
#                                      calendar='standard'))
#            else:
#                ts_x = f['PRECTOTLAND'][:]
#                ts_t = num2date(f['time'][:],
#                             units='hours since 1980-01-01 00:00',
#                             calendar='standard')
#        af.progress_bar(daf+1, len(da_file_list))
#    ts_x = ts_x*3600
#    return(lat, lon, ts_t, ts_x)


def get_da_ts(da_file):
    with Dataset(da_file, 'r') as f:
        lat = f['latitude'][:]
        lon = f['longitude'][:]
        ts_x = f['PRECTOTLAND'][:]
        ts_t = num2date(f['time'][:], units='hours since 1980-01-01 00:00',
                        calendar='standard')
        ts_x = ts_x*3600
    return({'lat': lat, 'lon': lon, 't': ts_t, 'x': ts_x})


def resampling_ts(input_ts, window, output_t):
    winsum = input_ts.rolling(window=window, min_periods=window).sum()
    resampled_ts = winsum[winsum.index.isin(output_t)]
    return(resampled_ts)
