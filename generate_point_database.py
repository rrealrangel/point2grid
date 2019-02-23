# -*- coding: utf-8 -*-
"""
Created on Thu Nov 02 21:13:36 2017

@author: r.realrangel
"""

# Load modules
from loc2utc import loc2utc
from read_vector_map import read_station_map
import auxfunx as af
import ConfigParser
import datetime as dt
import get_data_sets as getds
import io
import numpy as np
import os
import pandas as pd
import pickle
import pytz


# Load the configuration file
with open("config.ini") as f:
    configini = f.read()
cfg = ConfigParser.RawConfigParser(allow_no_value=True)
cfg.readfp(io.BytesIO(configini))


def import_station_map(vmap_name=cfg.get('inputs', 'station_map_file_name')):
    """Import the attributes of climatologic stations."""
    vmap = os.path.abspath(os.path.join(
            os.path.dirname('__file__'),
            '..',
            'inputs',
            vmap_name))
    return(read_station_map(vmap))


def cs_filename(num):
    return(cfg.get('inputs', 'cs_db_dir') + num + '.csv')


def merge_df(left, right):
    return(pd.merge(
            left,
            right,
            how='outer',
            left_index=True,
            right_index=True))


if __name__ == '__main__':
    af.open_close_message('open')
    station_attributes = import_station_map()
    curr = 0

    for num, atts in sorted(station_attributes.items()):
        # Generate the daily time series
        cs_ts_d = loc2utc(getds.get_cs_ts(cs_filename(num)), atts)
        cs_ts_d = cs_ts_d[cs_ts_d.index >= dt.datetime(
                1980, 1, 1).replace(tzinfo=pytz.utc)]

        if curr == 0:
            merged = cs_ts_d

        else:
            merged = merge_df(merged, cs_ts_d)

        curr += 1
        af.progress_bar(curr, len(station_attributes))

    outFileFullPath = os.path.abspath(os.path.join(
            os.path.dirname('__file__'), '..', 'inputs', 'tseries_1980'))
#    # Remove stations without not enough data since 1940.
#    dropStations = [sta for sta in merged.keys()
#                    if (1 - np.sum(np.isfinite(merged[sta])) /
#                        float(len(merged[sta])) >=
#                        float(cfg.get('inputs', 'maxMissedData')))]
#    merged = merged.drop(dropStations, axis=1)

    # Export results
    with open(outFileFullPath + '.obj', 'w') as objFile:
        pickle.dump(merged, objFile)

    merged.to_csv(outFileFullPath + '.csv')
