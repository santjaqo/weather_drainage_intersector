# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 17:36:05 2014

@author: S. Gaitan
"""
#you need to be working in 'xband_sobek' directory
import pandas as pd
import os
import fiona
from shapely.geometry import Point, Polygon, LineString, shape, mapping
import numpy as np
from rtree import index 
import datetime
#%%
def rainfall_loader(pathname, separated=','):
    '''Returns dataframe of rainfall event. Cell series are in rows, timeseries in columns.
    
    Arguments:    
    
    :pathname: --csv file, provided by Susana.
    
    The index of the returned dataframe is the same as the 'PROFILE' field in the rainfall cells.

    '''
    rainfall = pd.read_csv(pathname, sep=separated)
    rainfall = rainfall.transpose()
    rainfall = rainfall.drop(rainfall.columns[[0,1,2,3,4]], axis=1)#this deletes metadata columns, leaving only the rainfall data. There are 3600 columns referring to each of the cells of the 60*60 grid; and 118 rows, correspoinding to 118 time-steps.
    rainfall.index = [i + 1 for i in range(len(rainfall.index))]
    for i in range(len(rainfall.columns)):
        rainfall.iloc[:,i] = rainfall.iloc[:,i].apply(np.float)
    rainfall = rainfall / 60.0
    return rainfall
#%%
def links_centroider(links, name_for_centroids):
    with fiona.open(links) as input_shp:
    # change only the geometry of the schema: LineString -> Point
        input_shp.schema['geometry'] = "Point"
        # write the Point shapefile
        with fiona.open(name_for_centroids, 'w', 'ESRI Shapefile', input_shp.schema.copy(), input_shp.crs) as output_shp:
           for elem in input_shp:
               # GeoJSON to shapely geometry
               geom = shape(elem['geometry'])
               # shapely centroid to GeoJSON
               elem['geometry'] = mapping(geom.centroid)
               output_shp.write(elem)
#%%
def shp_writer(model_shp, geoDF, output_shp):
    with fiona.open(model_shp) as source:
        source_driver = source.driver
        source_crs = source.crs
        source_schema = source.schema
        #previous fields in properties are deleted
        del source_schema['properties']
        #a new field is set with its respectie data type
        source_schema['properties'] = {'mxint15min': 'float'}
        #writing a new file    
        with fiona.open(output_shp,
                        'w',
                        driver=source_driver,
                        crs=source_crs,
                        schema=source_schema) as collection:
            #rec = {'geometry': mapping(geoDF.loc[0].polygon),'properties':{'mxrai15min': 0.5}}
            #collection.write(rec)
            for i in geoDF.index:
                #create a record
                rec = {}
                #fill geometry
                rec['geometry'] = mapping(geoDF.loc[i].polygon)
                #fill attribute values
                intensity = float(geoDF.loc[i].maxrain_15min)
                rec['properties'] = {'mxint15min': intensity}
                collection.write(rec)
#%%
def fishnet_loader(a_shp):
    '''Returns a DataFrame containing the WKT polygons of the passed fishnet, and its crs.
    This is, basically, KNMI's 1Km^2 radar grid.
    Arguments

    :a_shp: --string path, the name of the shp file to be read.

    '''
    collection = fiona.open(a_shp, 'r')
    crs = collection.crs
    records = []
    for i in range(len(collection)):
        records.append( next(collection))
    collection.close()
    geoDF = pd.DataFrame({'type': [i['geometry']['type'] for i in records],
                          'properties': [i['properties'] for i in records],
                          'coordinates': [i['geometry']['coordinates'] for i in records]
                          }, index = [i['id'] for i in records])
    geoDF['sequence'] = np.arange(len(geoDF))
    return geoDF, crs
#%%
def fishing_points(a_fishnet, nodes):
    '''Returns arrays containing shapely polygons representing the passed rainfall_dataframe.
    
    Arguments:    
    
    :a_fishnet: --pd.DataFrame, created with the fishnet function.
    
    :a_ranfall_dataframe: --pd.DataFrame, created with the rainfall_dataframe function.

    '''
    fishnet_index = index.Index()
    geo_polygons = a_fishnet.coordinates.apply(lambda x: Polygon(x[0])).values
    for position, cell in enumerate(geo_polygons[:]):
        fishnet_index.insert(position, cell.bounds)
    def intersecter(x):
        return list(fishnet_index.intersection(Point(x[0], x[1]).bounds))
    def get_content(x):
        if x != []:
            return x[0]
        else:
            return -1
    a_rainfall['fishnet_index'] = a_rainfall['coordXY'].apply(lambda x: intersecter(x)).apply(lambda x: get_content(x))
    fished = a_rainfall[a_rainfall.fishnet_index != -1].merge(a_fishnet, how='left', left_on='fishnet_index', right_on='sequence', suffixes=('','_grid'))
    coordXY = fished.coordinates.apply(lambda x: Polygon(x[0]).exterior.coords[:])
    fished = fished.iloc[:,:-5]
    fished['coordXY'] = coordXY
    return fished
#%%
def nodes_indexes_loader(pathname, separated=','):
    nodes = pd.read_csv(pathname)
    nodes = nodes.drop(nodes.columns[0], axis=1)
    return nodes
#%%
#def rainfall_loader(pathname, resolution, separated=','):
#    '''Returns dataframe of rainfall event. Cell series are in rows, timeseries in columns.
#    
#    Arguments:    
#    
#    :pathname: --csv file, provided by Susana.
#    
#    :resolution: --int, the ammount of gridcells, e.g. 3600 for the 1ha resolution rainfall dataset.
#    
#    The index of the returned dataframe is the same as the 'PROFILE' field in the rainfall cells.
#
#    '''
#    rainfall = pd.read_csv(pathname, sep=separated)
#    rainfall = rainfall.transpose()
#    rainfall = rainfall.drop(rainfall.columns[[0,1,2,3,4]], axis=1)#this deletes metadata columns, leaving only the rainfall data. There are 3600 columns referring to each of the cells of the 60*60 grid; and 118 rows, correspoinding to 118 time-steps.
#    rainfall.index = [i + 1 for i in range(resolution)]
#    return rainfall
#%%
def index_joiner(nodes, rainfall):
    joined = pd.merge(nodes, rainfall, how='left', left_on=nodes.columns[-1], right_index=True)
    return joined

#if __name__ == "__main__":
    saved_path = os.getcwd()
    os.chdir('xycoordinatesnodes')
    nodes = nodes_indexes_loader('nodes_20110628_Kralingen.csv')
    os.chdir(saved_path)
    os.chdir('05_Data_Partners_ToSend/20110628_Conv')
    rainfall = rainfall_loader('RainGain_SimS_1min_100m_2011-06-28.CSV')
    os.chdir(saved_path)
    joined = index_joiner(nodes, rainfall)
#    os.chdir(os.path.dirname(os.getcwd()))
#    os.chdir('data_Adam_MetropolitanSolutions/radar_adam_MS/RAD_NL25_RAC_MFBS_5min')
#    wrap = resampler(wrapper(time_series, a_fishnet))
#    os.chdir(saved_path)
#%%

    rainfall = rainfall_loader('RainGain_SimS_1min_100m_2011-06-28.CSV')
    collection = fiona.open('kralingen_n_20110628_100m.shp')
    records = list(collection)
    collection.close()
    stations = pd.DataFrame({'PROFILE':[i['properties']['PROFILE'] for i in records]}, index= [i['properties']['ID'] for i in records])
    start_time = datetime.datetime(2011, 06, 28, 22,01)
    time_series = pd.date_range(start_time, periods=len(rainfall.columns), freq='60S')
    end_time = time_series[-1].to_datetime()
    joined = pd.merge(stations, rainfall, how='left', left_on=stations.columns[-1], right_index=True)
    joined.columns = [joined.columns[0]] + [i.to_datetime() for i in time_series]
    #creating the final .BUI file
    name_output_file = '20110628_100m_kralingen.BUI\n'
    a_file = open('03_10_14_example_meteofile.BUI')
    lines = []
    for i in a_file:    lines.append(i)
    lines[0] = lines[0][:lines[0].find('\\')] + name_output_file
    lines[1] = lines[1][:lines[1].find(':') + 2] + datetime.datetime.now().strftime('%d-%m-%Y') + '   ' + datetime.datetime.now().strftime('%H:%M:%S') + '\n'
    lines[6] = ' {} \n'.format(len(stations.columns))
    rest_of_lines = lines[11:]
    del lines[8:]
    for i in stations.index:    lines.append('Station_' + str(i) + '\n')
    rest_of_lines[2] = ' 1  60 \n'
    delta_time = datetime.datetime(1,1,1) + datetime.timedelta(seconds = (end_time - start_time).total_seconds()) #beware: year, month, and day have an initial dummy value of 1. It must be substracted from the day value. Days are limited to 31.
    rest_of_lines[7] = ' ' + start_time.strftime('%Y %m %d %H %M %S') + ' ' + '{} {} {} {}'.format((delta_time.day -1), delta_time.hour, delta_time.minute, delta_time.second) + '\n'
    del rest_of_lines[8:]
    for i in rest_of_lines:    lines.append(i)
    #writing rainfall values
    stations_to_write = joined.iloc[:,1:].transpose()
    for i in range(len(stations_to_write.index)):
        step = stations_to_write.iloc[i].values
        step_write = ' '
        for value in step:
            step_write += ('  ' + str(value))
        step_write += '\n'
        lines.append(step_write)
    #dumping the lines into a .BUI file:
    output_file = open(name_output_file[:-1], mode = 'w')
    output_file.writelines(lines)
    output_file.close()
    
    
        
    
    
    
    
    
    
    with fiona.open('fishnet_depth.shp') as source:
        source_driver = source.driver
        source_crs = source.crs
        #source_schema = source.schema
        #previous fields in properties are deleted
        #del source_schema['properties']
        source_schema = {'geometry': 'Polygon',
                         'properties': {'a_fieldname': 'float'}}
        #writing a new file    
        with fiona.open('output_shp.shp','w',driver=source_driver,crs=source_crs,schema=source_schema) as collection:
            rec = {}
            rec['geometry'] = mapping(Polygon([(0,0),(1,0),(1,1),(0,0)]))
            rec['properties'] = {'a_fieldname': float(np.float(3.0))}
            collection.write(rec)
    #%%
    with fiona.open('fishnet_depth.shp') as source:
        source_driver = source.driver
        source_crs = source.crs
        #source_schema = source.schema
        #previous fields in properties are deleted
        #del source_schema['properties']
        source_schema = {'geometry': 'Polygon',
                         'properties': {'a_fieldname': 'float'}}
        #writing a new file    
        with fiona.open('output_shp.shp','w',driver=source_driver,crs=source_crs,schema=source_schema) as collection:
            rec = {}
            rec['geometry'] = mapping(Polygon([(0,0),(1,0),(1,1),(0,0)]))
            rec['properties'] = {'a_fieldname': float(np.float(3.0))}
            collection.write(rec)

    