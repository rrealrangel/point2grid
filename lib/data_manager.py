#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""point2grid
Converts a set of ground-based observations into a gridded dataset. It
has been designed to work with the raw files of climatological
observations from the Climatologic Database of the National
Meteorological Service of Mexico (SMN).

AUTHOR
    Roberto A. Real-Rangel (Institute of Engineering UNAM; Mexico)

LICENSE
    GNU General Public License
"""
__author__ = 'Roberto A. Real-Rangel (Institute of Engineering UNAM)'
__license__ = 'GNU General Public License version 3'

import numpy as np
import sys

from pathlib2 import Path
from pyproj import Proj, transform
import toml
import xarray as xr

import lib.quality_control_tests as qct


class Configurations():
    """
    """
    def __init__(self, config_file):
        self.config_file = config_file

        with open(self.config_file, 'rb') as file_content:
            config = toml.load(file_content)

        for key, value in config.items():
            setattr(self, key, value)


class GridDataset:
    def __init__(self, xmin, xmax, ymin, ymax, xres, yres):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.xres = xres
        self.yres = yres
        self.x = np.round(np.arange(
                self.xmin + (self.xres / 2), self.xmax, self.xres), 4)
        self.y = np.round(np.arange(
                self.ymin + (self.yres / 2), self.ymax, self.yres), 4)
        self.values = np.zeros((len(self.y), len(self.x))) * np.nan

    def resample(self, new_grid):
        output = np.zeros((len(new_grid.y), len(new_grid.x))) * np.nan
        for i, x in enumerate(new_grid.x):
            for j, y in enumerate(new_grid.y):
                output[j, i] = np.nanmean(
                        [self.values[m, n]
                            for m in np.where(
                                    (self.y >= (y - (new_grid.yres / 2))) &
                                    (self.y <= (y + (new_grid.yres / 2))))[0]
                            for n in np.where(
                                    (self.x >= (x - (new_grid.xres / 2))) &
                                    (self.x <= (x + (new_grid.xres / 2))))[0]])

        return(output)

    def filter_by_stations(self, stations_x, stations_y, min_stations):
        output = np.zeros((len(self.y), len(self.x))) * np.nan

        for i, x in enumerate(self.x):
            for j, y in enumerate(self.y):
                if np.sum((stations_x >= (x - (self.xres/2))) &
                          (stations_x <= (x + (self.xres/2))) &
                          (stations_y >= (y - (self.yres/2))) &
                          (stations_y <= (y + (self.yres/2)))) >= min_stations:
                    output[j, i] = self.values[j, i]

        return(output)


def list_files(parent_dir, ext):
    """List all files in a directory with a specified extension.

    Parameters
        parent_dir: string
            Full path of the directory of which the files are to be listed.
        ext: string or list of strings
            Extension(s) of the files to be listed.
    """
    parent_dir = Path(parent_dir)
    files_list = []

    if isinstance(ext, str):
        patterns = ['**/*' + ext]

    elif isinstance(ext, list):
        patterns = ['**/*' + i for i in ext]

    for patt in patterns:
        files_list.extend(parent_dir.glob(pattern=patt))

    return(files_list)


def get_elevation(catalog_dir, inregion_attrs):
    catalog_paths = list_files(parent_dir=catalog_dir, ext=['.nc', '.nc4'])

    for station in inregion_attrs.keys():
        for path in catalog_paths:
            if station in str(path):
                inregion_attrs[station]['Elevation'] = xr.open_dataset(
                        path).attrs['Elevation']

    return(inregion_attrs)


def progress_bar(current, total, message="- Processing"):
    progress = float(current)/total
    sys.stdout.write("\r    {} ({:.1f} % processed)".format(
            message, progress * 100))

    if progress < 1:
        sys.stdout.flush()

    else:
        sys.stdout.write('\n')


def gen_database(catalog_dir, inregion_attrs, variable, qc):
    """ Read precipitation from all station files and generate a
    pandas.Dataframe with them.
    """
    catalog_paths = list_files(parent_dir=catalog_dir, ext=['.nc', '.nc4'])

    for station in inregion_attrs.keys():
        for path in catalog_paths:
            if station in str(path):
                inregion_attrs[station]['catalog_path'] = path

    inregion_data = xr.Dataset(coords={'time': [np.datetime64('today', 'D')]})

    for st, station in enumerate(inregion_attrs.values()):
        # Get the station data.
        time_series = xr.open_dataset(station['catalog_path'])[variable]
        time_series = time_series.to_dataset()
        min_time = min([min(time_series.time.values), min(inregion_data.time.values)])
        max_time = max([max(time_series.time.values), max(inregion_data.time.values)])
        new_index = np.arange(min_time, max_time, np.timedelta64(1, 'D'))
        inregion_data = inregion_data.reindex(time=new_index)

        # Apply quality control tests.
        if qc['apply']:
            # Perform tests to raw data.
            if qc['gross_range_test']:
                time_series[variable + 'gross_range_test'] = qct.range_test(
                        input_ts=time_series[variable],
                        threshold=qc['std_thresh'],
                        climatology=qc['climatology'])

            if qc['climatology_test']:
                time_series[variable + '_climtest'] = qct.range_test(
                        input_ts=time_series[variable],
                        threshold=qc['std_thresh'],
                        climatology=qc['climatology'])

            if qc['spikes_data_test']:
                time_series[variable + '_spiktest'] = qct.spikes_data_test(
                        input_ts=time_series[variable],
                        threshold=qc['std_thresh'],
                        climatology=qc['climatology'])

            if qc['change_rate_test']:
                time_series[variable + '_chrttest'] = qct.change_rate_test(
                        input_ts=time_series[variable],
                        threshold=qc['std_thresh'],
                        climatology=qc['climatology'])

            if qc['flat_series_test']:
                time_series[variable + '_flattest'] = qct.flat_series_test(
                        input_ts=time_series[variable],
                        value_tolerance=0.0,
                        repetitions_tolerance=2,
                        skipzero=qc['climatology'])

            if qc['tmp_outlier_test']:
                time_series[variable + '_outltest'] = qct.tmp_outlier_test(
                        input_ts=time_series[variable],
                        c=7.5,
                        threshold=qc['std_thresh'])

            # Remove suspicious values.
            tests = [i for i in time_series.var().keys() if i != variable]
            is_outlier = time_series[tests[0]]

            for test in tests:
                is_outlier = (is_outlier | time_series[test])

            time_series[variable][is_outlier] = np.nan

        # Remove from database stations with less than 10 years of data.
        if float(time_series[variable].notnull().sum() / 365.) < 10:
            pass

        else:
            # TODO: Perform a test to filter too porous datasets. This
            # test will not take into account years without data.
            # gaps = float(dataset['precipitation'].isnull().sum(dim='time'))
            # (gaps / dataset['precipitation'].size) > threshold
            inregion_data[station['ID']] = time_series[variable].copy()
            progress_bar(
                    current=(st + 1), total=len(inregion_attrs),
                    message=("- Retrieving precipitation data from"
                             " climatological stations in the region"))

    return(inregion_data.to_dataframe())
