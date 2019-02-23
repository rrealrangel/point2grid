# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 14:39:54 2018
@author: RRealR
Based on this code: https://goo.gl/e7Kp4X
"""

# For more info about the geostatsmodels module, go to https://bit.ly/2vFNGOx
from netCDF4 import Dataset, date2num
from pathlib2 import Path
from pyproj import Proj, transform
from read_vector_map import read_station_map, read_elevation_map
from scipy.spatial import distance_matrix
import gc
import numpy as np
import pandas as pd
import sys
import time


# Set some initial parameters
station_map_file_name = Path(
        'C:/Users/rreal/OneDrive/documentos/proyectos/2017/'
        'evaluacion_gldas_merra/analisis/precipitation/inputs/'
        'climatological_stations.shp')
stations_data_file_name = Path(
        'C:/Users/rreal/OneDrive/documentos/proyectos/2017/'
        'evaluacion_gldas_merra/analisis/precipitation/inputs/'
        'precipitation_database_1980-2015.csv')
elevation_map_file_name = Path(
        'C:/Users/rreal/OneDrive/documentos/proyectos/2017/'
        'evaluacion_gldas_merra/analisis/precipitation/inputs/'
        'Republica30_R100md.shp')
maxMissedData = 0.20


class Grid:
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


def progress_bar(current, total, item_id=''):
    progress = float(current)/total

    if item_id == '':
        sys.stdout.write("\r        Processing item {} of {} ({:.1f} %)".
                         format(current, total, progress * 100))

    else:
        sys.stdout.write("\r        Processing {} (progress: {:.1f} %)".
                         format(item_id, progress * 100))

    if progress < 1:
        sys.stdout.flush()

    else:
        sys.stdout.write('\n')


def station_location(station_time_series, station_attributes, elevation):
    stations_x = np.array(
            [station_attributes[i]['LONGITUDE'] for i in
             station_time_series.keys()])
    stations_y = np.array(
            [station_attributes[i]['LATITUDE'] for i in
             station_time_series.keys()])
    stations_z = np.array(
            [station_attributes[i]['ELEVATION'] for i in
             station_time_series.keys()])
    return(stations_x[np.isfinite(station_cross_section)],
           stations_y[np.isfinite(station_cross_section)],
           stations_z[np.isfinite(station_cross_section)])


def interpolate_idw(near_stations, station_value, power=2):
    nominator = 0
    denominator = 0

    for i in range(len(near_stations)):
        value = station_value.loc[near_stations.index[i]]
        distance = near_stations[i]
        nominator += (value / pow(distance, power))
        denominator += (1 / pow(distance, power))

    # Return NODATA if the denominator is zero
    if denominator > 0:
        return(nominator/denominator)
    else:
        return(np.nan)


def create_nc4(output_file, input_time, input_lat, input_lon, input_values):
    with Dataset(output_file, 'w', format='NETCDF4') as rootgrp:
        # rootgrp = Dataset(out, 'w')
        rootgrp.History = 'Created on ' + time.ctime(time.time())
        rootgrp.Source = 'Institute of Engineering UNAM'
        rootgrp.Description = ("Ground-based precipitation")
        rootgrp.createDimension('time', len(input_time))
        rootgrp.createDimension('lat', len(input_lat))
        rootgrp.createDimension('lon', len(input_lon))

        # Define time variable
        t = rootgrp.createVariable(
                varname='time', datatype='f4', dimensions=('time',), zlib=True,
                fill_value=-32768)
        t.units = 'days since 1980-01-01'
        t[:] = input_time

        # Define latitude (Y coordinate) variable
        lat = rootgrp.createVariable(
                varname='lat', datatype='f4', dimensions=('lat',), zlib=True,
                fill_value=np.nan)
        lat.units = 'Degrees north'
        lat[:] = input_lat

        # Define longitude (X coordinate) variable
        lon = rootgrp.createVariable(
                varname='lon', datatype='f4', dimensions=('lon',), zlib=True,
                fill_value=np.nan)
        lon.units = 'Degrees east'
        lon[:] = input_lon

        # Define main variable
        precipitation = rootgrp.createVariable(
                varname='precipitation', datatype='f4',
                dimensions=('time', 'lat', 'lon',), zlib=True,
                fill_value=np.nan)
        precipitation.units = 'mm'
        precipitation[:] = input_values


if __name__ == "__main__":
    # Import initial data.
    station_time_series = pd.read_csv(
        filepath_or_buffer=str(stations_data_file_name), index_col=0,
        parse_dates=True)
    station_attributes = read_station_map(str(station_map_file_name))
    elevation = read_elevation_map(str(elevation_map_file_name))

    # Reproject coordinates to get distances in meters (not in degrees).
    elevation['centroid_x_lcc'] = np.empty(
            shape=(np.shape(elevation['centroid_x']))) * np.nan
    elevation['centroid_y_lcc'] = np.empty(
            shape=(np.shape(elevation['centroid_y']))) * np.nan

    for i in range(len(elevation['elevation'])):
        elevation['centroid_x_lcc'][i], elevation['centroid_y_lcc'][i] = (
                transform(Proj(init='EPSG:4326'), Proj(init='EPSG:6362'),
                          elevation['centroid_x'][i],
                          elevation['centroid_y'][i]))

    date_groups = [(year, month)
                   for year
                   in list(set([date.year
                                for date
                                in station_time_series.index]))
                   for month
                   in range(1, 13)]

    for group in date_groups:
        year = group[0]
        month = group[1]
        print("Interpolating records of {}-{}".format(
                year, str(month).zfill(2)))
        auxiliar_base = []
        auxiliar_gldas25 = []
        auxiliar_gldas10 = []
        auxiliar_merra = []
        group_indices = station_time_series.loc[
                    str(year)+'-'+str(month)].index

        for count, date in enumerate(group_indices):
            progress_bar(count, len(group_indices))

            # Creating the output grid and populating the output matrix value
            base = Grid(-118, -86, 14, 34, 0.1, 0.1)
            gldas25 = Grid(-118, -86, 14, 34, 0.25, 0.25)
            gldas10 = Grid(-118, -86, 14, 34, 1.0, 1.0)
            merra = Grid(-118.4375, -86, 13.75, 34, 0.625, 0.5)

            # Get available data for the given date.
            station_cross_section = station_time_series.xs(date)
            stations_x, stations_y, stations_z = station_location(
                    station_time_series, station_attributes, elevation)
            stations_x_lcc = np.zeros(stations_x.shape) * np.nan
            stations_y_lcc = np.zeros(stations_y.shape) * np.nan
            for i in range(len(stations_x)):
                stations_x_lcc[i], stations_y_lcc[i] = transform(
                        Proj(init='EPSG:4326'), Proj(init='EPSG:6362'),
                        stations_x[i], stations_y[i])
            station_cross_section = station_cross_section[
                    np.isfinite(station_cross_section)]

            # Generate a distance matrix for the grid and the data.
            stations_points = zip(
                    stations_x_lcc,
                    stations_y_lcc,
                    stations_z)
            grid_points = zip(
                    elevation['centroid_x_lcc'],
                    elevation['centroid_y_lcc'],
                    elevation['elevation'])
            distmat = pd.DataFrame(
                    distance_matrix(x=grid_points, y=stations_points),
                    columns=station_cross_section.index)
            max_distance = 50000   # 50 km

            # Apply the IDW interpolation.
            for i in range(len(distmat)):
                x = elevation['centroid_x'][i]
                y = elevation['centroid_y'][i]
                base_x_index = np.where(base.x == x)[0]
                base_y_index = np.where(base.y == y)[0]
                cell_distances = distmat.iloc[i]
                near_stations = cell_distances[(cell_distances < max_distance)]
                base.values[base_y_index, base_x_index] = interpolate_idw(
                        near_stations, station_cross_section, power=2)

            gldas25.values = base.resample(new_grid=gldas25)
            gldas25.values = gldas25.filter_by_stations(
                    stations_x, stations_y, min_stations=1)
            gldas10.values = base.resample(new_grid=gldas10)
            gldas10.values = gldas10.filter_by_stations(
                    stations_x, stations_y, min_stations=1)
            merra.values = base.resample(new_grid=merra)
            merra.values = merra.filter_by_stations(
                    stations_x, stations_y, min_stations=1)

            auxiliar_base.append(np.expand_dims(base.values, axis=0))
            auxiliar_gldas25.append(np.expand_dims(gldas25.values, axis=0))
            auxiliar_gldas10.append(np.expand_dims(gldas10.values, axis=0))
            auxiliar_merra.append(np.expand_dims(merra.values, axis=0))
            progress_bar(count+1, len(group_indices))

        input_time = date2num(
                [i.replace(tzinfo=None)
                 for i in group_indices.to_pydatetime()],
                'days since 1980-01-01')

        # Export base datasets.
        output_dir = Path('base') / str(year)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = (output_dir / ('GBO.Pre.010.' + str(year) +
                                     str(month).zfill(2) + '.nc4'))
        create_nc4(output_file, input_time, base.y, base.x,
                   np.vstack(auxiliar_base))

        # Export GLDAS25 datasets.
        output_dir = Path('gldas25') / str(year)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = (output_dir / ('GBO.Pre.G225.' + str(year) +
                                     str(month).zfill(2) + '.nc4'))
        create_nc4(output_file, input_time, gldas25.y, gldas25.x,
                   np.vstack(auxiliar_gldas25))

        # Export GLDAS10 datasets.
        output_dir = Path('gldas10') / str(year)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = (output_dir / ('GBO.Pre.G210.' + str(year) +
                                     str(month).zfill(2) + '.nc4'))
        create_nc4(output_file, input_time, gldas10.y, gldas10.x,
                   np.vstack(auxiliar_gldas10))

        # Export MERRA2 datasets.
        output_dir = Path('merra') / str(year)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = (output_dir / ('GBO.Pre.M2.' + str(year) +
                                     str(month).zfill(2) + '.nc4'))
        create_nc4(output_file, input_time, merra.y, merra.x,
                   np.vstack(auxiliar_merra))

        del(distmat, base, gldas25, gldas10, merra)
        gc.collect()   # Trying to free some memory
