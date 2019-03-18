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
import osr

import gdal
import ogr

from data_manager import GridDataset


def get_shapefile_fields(map_file):
    """Read the position of a climatological station from a vector map.

    Parameters
        map_file: string
            Full path of the input vector map file.

    Source
        Python GDAL/OGR Cookbook : Get Shapefile Fields - Get the user
        defined fields (https://bit.ly/2u8ZaXg)
    """
    data_source = ogr.Open(map_file)
    first_layer = data_source.GetLayer(0)
    layer_definition = first_layer.GetLayerDefn()
    fields = [layer_definition.GetFieldDefn(i).GetName()
              for i in range(layer_definition.GetFieldCount())]
    table = {}

    for feature in first_layer:
        table[feature.GetField('ID')] = {
                field: feature.GetField(field) for field in fields}

    return(table)


def create_buffer(input_vmap, buffer_dist):
    """
    Buffers features of a layer and saves them to a new Layer.
    Source: https://pcjericks.github.io/py-gdalogr-cookbook
    """
    # Open an input datasource
    input_source = ogr.GetDriverByName('ESRI Shapefile').Open(input_vmap)
    input_layer = input_source.GetLayer(0)

    # Create an output datasource in memory
    output_source = ogr.GetDriverByName('Memory').CreateDataSource('out')

    # Open the memory datasource with write access
    proj = osr.SpatialReference()
    proj.SetWellKnownGeogCS("EPSG:6372")
    output_layer = output_source.CreateLayer(
            'Buffer', geom_type=ogr.wkbPolygon, srs=proj)
    feature_definition = output_layer.GetLayerDefn()

    for feature in input_layer:
        input_geometry = feature.GetGeometryRef()
        output_geometry = input_geometry.Buffer(buffer_dist)
        output_feature = ogr.Feature(feature_definition)
        output_feature.SetGeometry(output_geometry)
        output_layer.CreateFeature(output_feature)
        output_feature = None

    return(output_source)


def spatial_filter(features, ref):
    """Return a list of rain gauges near to the point of interest.

    Parameters
        features: string
            Ful path of the vector map of interest.
        ref: osgeo.ogr.DataSource
            Reference region.

    Source
        Python GDAL/OGR Cookbook : Spatial Filter (https://bit.ly/
        2F0LhiU)
    """
    # Get (first) layer from the reference region.
    ref_feature = ref.GetLayer(0)[0]

    # Open interest features.
    features_source = ogr.GetDriverByName('ESRI Shapefile').Open(features)
    features_layer = features_source.GetLayer(0)

    # Create an output datasource in memory
    output_source = ogr.GetDriverByName('Memory').CreateDataSource('out')

    # Open the memory datasource with write access
    ogr.GetDriverByName('Memory').Open('out', 1)

    # Copy a layer to memory
    output_source.CopyLayer(features_layer, 'Stations', ['OVERWRITE=YES'])

    # Apply spatial filter.
    output_layer = output_source.GetLayer('Stations')
    output_layer.SetSpatialFilter(ref_feature.GetGeometryRef())

    return(output_source)


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


def inregion_cells(study_region, res, nodata):
    """
    Source:
        https://bit.ly/2HxeOng
    """
    # Open the data source and read in the extent
    source_ds = ogr.Open(study_region)
    source_layer = source_ds.GetLayer(0)
    xmin, xmax, ymin, ymax = source_layer.GetExtent()

    def round_mult(num, mult, to):
        if to == 'up':
            return(int(mult * round(float(num + mult) / mult)))

        elif to == 'down':
            return(int(mult * round(float(num - mult) / mult)))

    xmin = round_mult(num=xmin, mult=res, to='down')
    xmax = round_mult(num=xmax, mult=res, to='up')
    ymin = round_mult(num=ymin, mult=res, to='down')
    ymax = round_mult(num=ymax, mult=res, to='up')
    x_coords = np.arange((xmin + (res / 2)), xmax, res)
    y_coords = np.arange((ymin + (res / 2)), ymax, res)

    # Create the destination data source
    cols = int((xmax - xmin) / res)
    rows = int((ymax - ymin) / res)
    output_source = gdal.GetDriverByName('MEM').Create(
            '', cols, rows, gdal.GDT_Byte)
    output_source.SetGeoTransform((xmin, res, 0, ymax, 0, -res))
    output_band = output_source.GetRasterBand(1)
    output_band.SetNoDataValue(nodata)

    # Rasterize
    gdal.RasterizeLayer(output_source, [1], source_layer, burn_values=[1])

    # Read as array
    array = output_band.ReadAsArray().astype(int)
    cells = []

    for col, x_coord in enumerate(x_coords):
        for row, y_coord in enumerate(y_coords[::-1]):
            if array[row, col] == 1:
                cells.append([x_coord, y_coord])

    grid = GridDataset(
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, xres=res, yres=res)

    return(cells, grid)


def retrieve_elevation(input_rmap, coordinates, res, nodata):
    """
    Source:
        https://bit.ly/2XZe1kJ
    """
#    driver = gdal.GetDriverByName('GTiff')
#    filename = "/home/zeito/pyqgis_data/aleatorio.tif" #path to raster
    input_source = gdal.Open(input_rmap)
    input_band = input_source.GetRasterBand(1)
    cols = input_source.RasterXSize
    rows = input_source.RasterYSize
    (xmin, x_size, x_rotation, ymax, y_rotation, y_size) = (
            input_source.GetGeoTransform())
    y_size = y_size * -1
    data = input_band.ReadAsArray(0, 0, cols, rows).astype(float)
    data[data == nodata] = np.nan
    x_coords = np.arange(xmin, (xmin + (cols * x_size)), x_size)
    y_coords = np.arange((ymax - ((rows - 1) * y_size)), ymax + y_size, y_size)
    y_coords = y_coords[::-1]

    for p, point in enumerate(coordinates):
        irows = np.where((y_coords > (point[1] - (res / 2))) &
                         (y_coords < (point[1] + (res / 2))))[0]
        icols = np.where((x_coords > (point[0] - (res / 2))) &
                         (x_coords < (point[0] + (res / 2))))[0]
        coordinates[p].append(
                np.mean([data[i, j] for i in irows for j in icols]))

    return(coordinates)


def filter_stations(catalog_vmap, study_region, buffer_dist=80000):
    """Return a subset of climatological stations filtered by location.

    Filter a catalog of climatological stations of the Climatological
    Database (BDC) of the National Meteorological Service of Mexico
    (SMN) within an influen region defined by a buffer distance from a
    study region.

    Parameters:
        catalog_vmap : string
            Full path of the vector map (shapefile) that contains the
            catalog of climatolog stations of the BDC.
        study_region : string
            Full path of the vector map (shapefile) that contains the
            region of interest.
        buffer_dist : float or integer, optional
            The distance that reach the region of influence of the
            values observed in a climatological stations. By default,
            buffer_dist is 80000 (m; 80 km).

    Returns:
        inregion_vmap, inregion_attrs : osgeo.ogr.DataSource, dictionary
            inregion_vmap contains all the features (climatological
            stations) retained by the location filter. inregion_attrs
            contains the attributes of those features.
    """
    # List all stations available and their locations (lat, lon).
    catalog_attrs = get_shapefile_fields(map_file=catalog_vmap)

    # Get region of influence (buffer) of the catchment.
    buffer_region = create_buffer(
            input_vmap=study_region,
            buffer_dist=buffer_dist)

    # Filter datasets within the buffer region.
    inregion_vmap = spatial_filter(features=catalog_vmap, ref=buffer_region)
    inregion_attrs = {
            st.GetField('ID'): catalog_attrs[st.GetField('ID')]
            for st in inregion_vmap[0]}

    return(inregion_vmap, inregion_attrs)
