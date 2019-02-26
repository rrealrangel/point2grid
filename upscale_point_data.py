# -*- coding: utf-8 -*-
""" Rasteriza precipitation
    Author
    ------
        Roberto A. Real-Rangel (Institute of Engineering UNAM; Mexico)

    License
    -------
        GNU General Public License version 3
"""
from pathlib2 import Path
from pyproj import Proj, transform
from scipy.spatial import distance_matrix
import gdal
import numpy as np
import pandas as pd
import sys
import toml
import xarray as xr

with open('config.toml', 'rb') as config_file:
    config = toml.load(config_file)


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


def progress_bar(current, total, message="- Processing"):
    progress = float(current)/total
    sys.stdout.write("\r    {} ({:.1f} % processed)".format(
            message, progress * 100))

    if progress < 1:
        sys.stdout.flush()

    else:
        sys.stdout.write('\n')


def gen_database(input_dir):
    """ Read precipitation from all station files and generate a
    pandas.Dataframe whit them.
    """
    input_list = sorted(list(Path(input_dir).glob(pattern='**/*.nc')))
    prec_database = xr.Dataset()

    for inpf, input_file in enumerate(input_list):
        # Get the station data.
        vars2drop = [i for i in xr.open_dataset(input_file).var().keys()
                     if 'precipitation'not in i]
        dataset = xr.open_dataset(input_file, drop_variables=vars2drop)

        # Remove other variables not precipitation.
        tests = [i for i in dataset.var().keys() if i != 'precipitation']
        is_outlier = dataset[tests[0]]

        # Remove suspicious values.
        for test in tests:
            is_outlier = (is_outlier | dataset[test])

        dataset['precipitation'][is_outlier] = np.nan

        # Remove from database stations with less than 10 years of data.
        if float(dataset['precipitation'].notnull().sum() / 365.) < 10:
            pass

        else:
            # TODO: Perform a test to filter too porous datasets. This
            # test will not take into account years without data.
            # gaps = float(dataset['precipitation'].isnull().sum(dim='time'))
            # (gaps / dataset['precipitation'].size) > threshold

            prec_database[dataset.attrs['ID']] = (
                    dataset['precipitation'].copy())
            is_outlier = is_outlier.reindex(time=sorted(
                    prec_database.time.values)).astype('bool')
            progress_bar(
                    current=(inpf + 1), total=len(input_list),
                    message=("- Retrieving precipitation data from all "
                             "climatological stations"))

    return(prec_database.to_dataframe())


def interpolate_idw(distances, values, power=2):
    nominator = 0
    denominator = 0

    for i in range(len(distances)):
        value = values.loc[distances.index[i]]
        distance = distances[i]
        nominator += (value / pow(distance, power))
        denominator += (1 / pow(distance, power))

    # Return NODATA if the denominator is zero
    if denominator > 0:
        return(nominator/denominator)
    else:
        return(np.nan)


def retrieve_dem_elev(input_file, points_list, nodata=-32768):
    """ Retrieve pixel value with coordinates as input.
    Parameters
        input_file: string
        points_list: list of tuples

    Reference
        https://bit.ly/2VkQvML
    """
    # TODO: Retrieve the elevation of climatological stations from the
    # DEM, instead of the metadata.
    driver = gdal.GetDriverByName('GTiff')
    raster = gdal.Open(input_file)
    band = raster.GetRasterBand(1)
    cols = raster.RasterXSize
    rows = raster.RasterYSize
    (x_min, x_size, x_rotation, y_max, y_rotation, y_size) = (
            raster.GetGeoTransform())
    x_max = x_min + (cols * x_size)
    y_min = y_max - (rows * abs(y_size))
    data = band.ReadAsArray(0, 0, cols, rows)
    output = []

    for point in points_list:
        x = point[0]
        y = point[1]

        if (x_min < x < x_max) & (y_min < y < y_max):
            col = int((point[0] - x_min) / x_size)
            row = int((y_max - point[1]) / abs(y_size))

            if float(data[row][col]) != nodata:
                output.append([point[0], point[1], float(data[row][col])])

    return(output)


if __name__ == "__main__":
    station_time_series = gen_database(input_dir=config['inputs']['datasets'])
    attrs = {}
    input_list = sorted(list(Path(
            config['inputs']['datasets']).glob(pattern='**/*.nc')))

    for st, station in enumerate(input_list):
        code = xr.open_dataset(station).attrs['ID']
        attrs[code] = dict(xr.open_dataset(station).attrs)
        attrs[code]['X_epsg4326'] = attrs[code].pop('LONGITUDE')
        attrs[code]['X_epsg4326'] = float(attrs[code]['X_epsg4326'])
        attrs[code]['Y_epsg4326'] = attrs[code].pop('LATITUDE')
        attrs[code]['Y_epsg4326'] = float(attrs[code]['Y_epsg4326'])
        attrs[code]['Z_epsg4326'] = attrs[code].pop('ELEVATION')
        attrs[code]['Z_epsg4326'] = float(attrs[code]['Z_epsg4326'])
        (attrs[code]['X_epsg6372'],
         attrs[code]['Y_epsg6372'],
         attrs[code]['Z_epsg6372']) = transform(
                p1=Proj(init='EPSG:4326'),
                p2=Proj(init='EPSG:6372'),
                x=attrs[code]['X_epsg4326'],
                y=attrs[code]['Y_epsg4326'],
                z=attrs[code]['Z_epsg4326'])
        progress_bar(
                st + 1, len(input_list),
                message="- Transforming projections of stations coordinates")

    years = range(
            station_time_series.index[0].year,
            station_time_series.index[-1].year + 1)

    for yy, year in enumerate(years):
        aux = []
        group_indices = station_time_series.loc[str(year)].index

        for count, date in enumerate(group_indices):
            basin_grid = GridDataset(
                    xmin=config['region']['xmin'],
                    xmax=config['region']['xmax'],
                    ymin=config['region']['ymin'],
                    ymax=config['region']['ymax'],
                    xres=config['region']['xres'],
                    yres=config['region']['yres'])

            # Get available data for the given date.
            values = station_time_series.xs(date).dropna()

            # Generate a distance matrix for the grid and the data.
            stations_points = [
                    (attrs[i]['X_epsg6372'], attrs[i]['Y_epsg6372'],
                     attrs[i]['Z_epsg6372']) for i in values.keys()]
            grid_cells = retrieve_dem_elev(
                    input_file=config['inputs']['dem'],
                    points_list=[(x, y)
                                 for y in basin_grid.y
                                 for x in basin_grid.x])
            distmat = pd.DataFrame(
                    distance_matrix(x=grid_cells, y=stations_points),
                    columns=values.index)

            # Apply the IDW interpolation.
            for c, cell in enumerate(grid_cells):
                x_index = np.where(basin_grid.x == cell[0])[0]
                y_index = np.where(basin_grid.y == cell[1])[0]
                cell_distances = distmat.iloc[c]
                distances = cell_distances[
                        (cell_distances < config['analysis']['max_distance'])]

                if len(distances) > 0:
                    basin_grid.values[y_index, x_index] = interpolate_idw(
                            distances=distances,
                            values=values,
                            power=2)

            aux.append(np.expand_dims(basin_grid.values, axis=0))

        output_dataarray = xr.DataArray(
                data=np.vstack(aux),
                coords={'time': group_indices,
                        'north': basin_grid.y,
                        'east': basin_grid.x},
                dims=['time', 'north', 'east'])

        # Export basin_grid datasets.
        output_file = (config['outputs']['output_dir'] + '/sg30057_prec100m_' +
                       str(year) + '.nc4')
        output_dataarray.to_netcdf(
                path=output_file,
                mode='w',
                format='NETCDF4')
        progress_bar(current=(yy + 1), total=len(years), message=(
                "- Interpolating daily records and exporting annual datasets"))
