# -*- coding: utf-8 -*-
"""
Created on Sat Oct 25 13:24:39 2014

@author: S. Gaitan
"""
#%%
import pandas as pd
import os
import fiona
#from shapely.geometry import Point, Polygon, LineString, shape, mapping
import numpy as np
#from rtree import index 
import datetime
from shapely.geometry import Point, Polygon, shape
from rtree import index
#%%    
#Rainfall time series are processed as pd.DataFrames and pickled
def rainfall_loader(pathname, separated=','):
    '''Returns dataframe of rainfall event. Cell series are in rows, timeseries in columns.
    
    Arguments:    
    
    :pathname: --csv file, provided by Susana.
    
    The index of the returned dataframe is the same as the 'PROFILE' field in the rainfall cells.

    '''
    rainfall = pd.read_csv(pathname, sep=separated, skiprows=[0,1,2,3,4,5], index_col=0, parse_dates=True)
    del rainfall.index.name
    rainfall = rainfall.transpose()
    for i in range(len(rainfall.columns)):
        rainfall.iloc[:,i] = rainfall.iloc[:,i].apply(np.float)
    return rainfall
#%%
def rainfall_pickler():
    saved_path = os.getcwd()
    os.chdir('events')
    events_folder = os.getcwd()
    os.chdir(saved_path)
    os.chdir('CSV-per-event')
    folders = []
    for (dir_path, dir_names, file_names) in os.walk(os.getcwd()):
        folders += dir_names
        break
    for folder in folders:
        os.chdir(folder)
        storms = []
        for (dir_path, dir_names, file_names) in os.walk(os.getcwd()):
            storms += file_names
            #os.chdir(saved_path)
            #break
        for storm in storms:
            rainfall_loader(saved_path + '/CSV-per-event/' + folder + '/' + storm).to_pickle(events_folder + '/' + storm + '.pickled')
        os.chdir(saved_path + '/CSV-per-event')
#%%
def fishnet_loader(a_shp):
    '''Returns a list containing geojson-like the WKT polygons of the passed fishnet.
    Arguments

    :a_shp: --string path, the name of the shp file to be read.

    '''
    with fiona.open(a_shp) as collection:
        fishnet = list(collection.items())
    return fishnet
#%%
def intersecter(a_fishnet_index, geo_bounds):
        return list(a_fishnet_index.intersection(geo_bounds))
def get_content(a_intersect_list):
    if a_intersect_list != []:
        #One point is evaluated at a time, len of returned list is always 1.
        return a_intersect_list
    else:
        return -1
#%%
def bui_typer(some_rainfall_stations, rainfall, time_res, spati_res, event_name, direction):
    '''Writes the .bui file for a storm in the prefered storm direction.
    Arguments

    :some_rainfall_stations: --pd.DataFrame, dataframe of the stations to be written.
    :rainfall: --pd.DataFrame, rainfall loaded from pickled file
    
    '''
    start_time = rainfall.columns[0].to_datetime()
    end_time = rainfall.columns[-1].to_datetime()
    a_timestep = rainfall.columns[-1].to_datetime() - rainfall.columns[-2].to_datetime()
    #seconds per time_step, and number of timesteps to be extended at the end of the rainfall event, are defined:
    if time_res == '1min':
        seconds_per_timestep = 60
        extender_module = 240
    elif time_res == '3min':
        seconds_per_timestep = 180
        extender_module = 80
    elif time_res == '5min':
        seconds_per_timestep = 300
        extender_module = 48
    elif time_res == '10min':
        seconds_per_timestep = 600
        extender_module = 24
    # four hours of zero rainfall are added at the end of the event. This allows to check how the sewer network resumes after the rainfall.
    time_offset = pd.date_range(end_time + a_timestep, end_time + (a_timestep * extender_module), freq= str(seconds_per_timestep) + 'S')
    for i in time_offset:
        some_rainfall_stations[i] = np.zeros(len(some_rainfall_stations.index))
    #end_time is updated:
    end_time = some_rainfall_stations.columns[-1].to_datetime()
    #a time delt with sobek format is written:
    delta_time = datetime.datetime(1,1,1) + datetime.timedelta(seconds = (end_time - start_time).total_seconds()) #beware: year, month, and day have an initial dummy value of 1. It must be substracted from the day value. Days are limited to 31.
    name_output_file = event_name + '_' + spati_res + '_' + time_res + '_' + direction + '.BUI'
    a_file = open('example_meteofile.BUI')
    lines = []
    for i in a_file:    lines.append(i)
    lines[0] = lines[0][:lines[0].find('\\')] + name_output_file + '\n'
    lines[1] = lines[1][:lines[1].find(':') + 2] + datetime.datetime.now().strftime('%d-%m-%Y') + '   ' + datetime.datetime.now().strftime('%H:%M:%S') + '\n'
    lines[6] = ' {} \n'.format(len(some_rainfall_stations))
    rest_of_lines = lines[11:]
    del lines[8:]
    #names of rainfall stations are written:
    for i in some_rainfall_stations.index:    lines.append("'Station_" + str(i) + "'\n")
    #seconds per timestep are written:
    rest_of_lines[2] = ' 1  {} \n'.format(seconds_per_timestep)
    #Eerste record bevat startdatum en -tijd, lengte van de gebeurtenis in dd hh mm ss. Het format is: yyyymmdd:hhmmss:ddhhmmss:
    rest_of_lines[7] = ' ' + start_time.strftime('%Y %m %d %H %M %S') + ' ' + '{} {} {} {}'.format((delta_time.day -1), delta_time.hour, delta_time.minute, delta_time.second) + '\n'
    del rest_of_lines[8:]
    for i in rest_of_lines:    lines.append(i)
    stations_to_write = some_rainfall_stations.iloc[:,1:].transpose()
    #rainfall series are written:    
    for i in range(len(stations_to_write.index)):
        step = stations_to_write.iloc[i].values
        step_write = ' '
        for value in step:
            step_write += ('  ' + str(value))
        step_write += '\n'
        lines.append(step_write)
    #dumping the lines into a .BUI file:
    output_file = open(name_output_file, mode = 'w')
    output_file.writelines(lines)
    output_file.close()
#%%
def joiner(direction):
    '''Joins the nodes with the pickled events and respective grid in the given direction. Must be run in the folder containing the events, grid, and node directories.
    Arguments

    :direction: --string, either 'perpendicular' or 'paralel'

    '''
    saved_path = os.getcwd()    
    os.chdir('pluvius_files/example')#Pluvius file is loaded. It is used to select only the components receiving surface water.
    a_file = open('PLUVIUS.3B')
    lines = [] #registries of pluvius file are stored in a list
    for i in a_file:
        lines.append(i)
    a_file.close()
    runoff_nodes = [] #the names of the components are stored in a separated list
    for i in lines:
        ini_slice = i.find("'")
        end_slice = i.find("'", ini_slice + 1)
        runoff_nodes.append(i[ini_slice + 1 : end_slice])    
    os.chdir(saved_path)
    os.chdir('nodes')#entering nodes folder. 
    with fiona.open('kralingen_n.shp') as collection:
        raw_nodes = list(collection.items())
    # only the nodes receiving rainfall runoff are consideraded
    nodes = []
    for node in raw_nodes:
        if node[1]['properties']['ID        '] in runoff_nodes:
            nodes.append(node)
    os.chdir(saved_path)
    os.chdir('events')#entering events folder.
    pickles = []#storing events into a list
    for (dir_path, dir_names, file_names) in os.walk(os.getcwd()):
        pickles += file_names
        break
    os.chdir(saved_path)
    #information for calling respective grids is extracted from the pickles name:
    #os.chdir(os.path.dirname(os.getcwd()))
    for storm in pickles:
        os.chdir(saved_path + '/events')
        rainfall = pd.read_pickle(storm)
        info_line = storm[(storm.find('S_')+2):(storm.find('.CSV'))]
        time_res, spati_res, event_name = info_line.split('_')
        #rainfall units are converted from mm/h to the adequate depth given the time_res
        #this is required by sobek
        if time_res == '1min':
            conversion_module = 60
        elif time_res == '3min':
            conversion_module = 20
        elif time_res == '5min':
            conversion_module = 12
        elif time_res == '10min':
            conversion_module = 6        
        rainfall = rainfall / conversion_module
        a_grid = 'grid_' + spati_res
        os.chdir(saved_path)        
        os.chdir('grids')
        #the respective rainfall grid is loaded:
        grid = fishnet_loader(a_grid + '_centrum' + '.shp')
        fishnet_index = index.Index()
        for cell in grid:        
            fishnet_index.insert(int(cell[1]['properties']['PROFILE']), shape(cell[1]['geometry']).bounds)
        #the fishnet_index is queried against the nodes
        nodes_indexed = []
        for node in nodes:
            geo_bounds = shape(node[1]['geometry']).bounds
            ###
            a_intersect_list = get_content(intersecter(fishnet_index, geo_bounds))
            for element in a_intersect_list:
                indexed_node = -1
                for cell in grid:
                    if element == cell[1]['properties']['PROFILE']: 
                        if shape(node[1]['geometry']).intersects(shape(cell[1]['geometry'])):
                            indexed_node = element
                    if indexed_node != -1:
                        break
                if indexed_node != -1:
                    break
            #node id and respective rainfall index are stored in a list:            
            node_id = node[1]['properties']['ID        ']
            nodes_indexed += [(node_id, str(indexed_node))]
        nodes_indexed = pd.DataFrame({'cell_id':[i[1] for i in nodes_indexed]}, index=[i[0] for i in nodes_indexed])
        #rainfall series are annotated to each node
        rainfall_stations = pd.merge(nodes_indexed, rainfall, how='left', left_on='cell_id', right_index=True)
        rainfall_stations = rainfall_stations[rainfall_stations.cell_id != -1]
        # check whether only run-off receiving nodes must be included in the .BUI
        # .bui file is written:
        os.chdir(saved_path)
        os.chdir('bui_files')
        bui_typer(rainfall_stations, rainfall, time_res, spati_res, event_name, direction)
        os.chdir(saved_path)
#%%
def pluvius_writer():
    saved_path = os.getcwd()    
    os.chdir('pluvius_files/example')#Pluvius file is loaded. It is used to select only the components receiving surface water.
    a_file = open('PLUVIUS.3B')
    lines = [] #registries of pluvius file are stored in a list
    for i in a_file:
        lines.append(i)
    a_file.close()
    new_lines = [] #the names of the components are stored in a separated list
    for i in lines:
        ini_slice = i.find("'")
        second_slice = i.find("'", ini_slice + 1)
        third_slice = i.find("'", second_slice + 1)
        fourth_slice = i.find("'", third_slice +1)
        new_lines.append(i[:third_slice +1] + 'Station_' + i[ini_slice + 1 : second_slice] + i[fourth_slice:])  
    os.chdir(saved_path)
    a_file = open('pluvius_files/example/PLUVIUS.3B', 'w')
    for i in new_lines:
        a_file.write(i)
    a_file.close()
    os.chdir(saved_path)
#%%
if __name__ == "__main__":
    #must be run in the folder containing the csvs, grids, nodes, and events.
    saved_path = os.getcwd()    
    #rainfall_pickler()
    os.chdir(saved_path)
    joiner('perpendicular')
    os.chdir(saved_path)
    #joiner('parallel') #Storm direction is not important in Rotterdam conditions
    pluvius_writer()