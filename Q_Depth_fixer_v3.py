# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 10:25:56 2015
Comma fixer: this code fixes the separation value of result files produced in Sobek simulations.
@author: S. Gaitan
"""
#%%
# THIS PROGRAM MUST BE RUN WITHIN THE SAME FOLDER WHERE THE SOBEK RESULT FILES ARE!
import os
import pandas as pd
import datetime
import multiprocessing as mp
import time
#%%
def reader(a_file):
    opened_file = open(str(a_file))
    lines = []
    counter = 0    
    for line in opened_file:
        counter += 1
        if counter == 5:
            header = line.strip('\n')
            header = header.split(',')
            header = header[1:]
        elif counter <= 4:
            continue
        else:
            lines.append(line)
    opened_file.close()
    return lines, header
#%%
def fixer(lines):
    new_lines = range(len(lines))
    for l_position, line in enumerate(lines):
        updated_line = []
        hit = 0
        length = len(line)
        for position, character in enumerate(line):
            if position == 0:
                continue
            elif position == length:
                continue
            elif character == ',':
                if line[position + 1] == ' ' or line[position + 1] == '-':
                    updated_line.append(line[hit:position])
                    hit = position + 1
        last_comma = line.rfind(' ,')
        updated_line.append(line[last_comma + 2 :].strip('\n'))
        updated_line = [value.replace(',', '.') for value in updated_line]
        updated_line[1:] = [float(value) for value in updated_line[1:]]
        new_lines[l_position] = updated_line
    return new_lines
#%%
def framer(new_lines, header):
    df = pd.DataFrame(new_lines)
    df = df.set_index(df.columns[0])
    df.columns = header
    df.index = [datetime.datetime.strptime(value, '%Y-%m-%d %H:%M:%S') for value in df.index]
    df.index.name = 'storm_time'
    return df
#%%
def writer(df, a_file):
    df.to_csv(a_file[:-4] + '_fxd.csv')
    #df.to_pickle(a_file[:-4] + '_fxd.pickled') #this writes the DataFrame object to a pickle. It is useful if you want to use the dataframe for analysis within Python.
#%%
def wrapper(a_file):
    lines, header = reader(a_file)
    new_lines = fixer(lines)
    df = framer(new_lines, header)
    writer(df, a_file)
#%%
if __name__ == "__main__":
    start_time = time.time()
    saved_path = os.getcwd()
    rslt_list = list(os.walk(saved_path))
    rslt_list = rslt_list[0][2]
    new_rslt_list = []
    for i in rslt_list:
        if i[-3:] != '.py' and i[-4:] != '.pyc':
            new_rslt_list.append(i)
    #rslt_list = rslt_list[:20]
    processors = mp.cpu_count()
#    pool = mp.Pool(processors)
#    pool.map(wrapper, rslt_list)
#Attempt with a Queue approach
    queue = mp.Queue()
    for i in new_rslt_list:
        queue.put(i)
    processes = [mp.Process(target = wrapper, args=(queue,)) for core in range(processors)]
    for core in processes:
        core.start()
    for core in processes:
        core.join()
    print("--- %s seconds ---" % (time.time() - start_time))
#%% 

     
    
    