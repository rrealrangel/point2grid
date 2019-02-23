# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 08:02:33 2018

@author: RRealR
"""
import numpy as np
import ogr


def read_station_map(full_path):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    vector_map = driver.Open(full_path)
    vector_map_layer = vector_map.GetLayer()
    vector_map_features = vector_map_layer.GetFeature(0)
    vector_map_layer.ResetReading()
    station_attributes = {}

    for feature in vector_map_layer:
        station_attributes[str(feature.GetField('CODE')).zfill(5)] = {
                'LATITUDE': feature.GetField('LATITUDE'),
                'LONGITUDE': feature.GetField('LONGITUDE'),
                'ELEVATION': feature.GetField('ELEVATION'),
                'TIME_ZONE': feature.GetField('TIMEZONE')}

    return(station_attributes)


def read_elevation_map(full_path):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    vector_map = driver.Open(full_path)
    vector_map_layer = vector_map.GetLayer()
    vector_map_features = vector_map_layer.GetFeature(0)
    vector_map_layer.ResetReading()
    centroid_x = np.array([feature.GetField('lon') for feature in
                           vector_map_layer], dtype=float)
    vector_map_layer.ResetReading()
    centroid_y = np.array([feature.GetField('lat') for feature in
                           vector_map_layer], dtype=float)
    vector_map_layer.ResetReading()
    elevation = np.array([feature.GetField('elevation') for feature in
                          vector_map_layer], dtype=float)
    return({'centroid_x': centroid_x, 'centroid_y': centroid_y,
            'elevation': elevation})
