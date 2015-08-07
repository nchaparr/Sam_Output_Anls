import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from matplotlib import rcParams
rcParams.update({'font.size': 10})
import pandas as pd
import glob
from collections import defaultdict
import re
import os
import h5py

#split 'theta_bar0000003420' into two parts
expr=re.compile('(.*)?(\d{10,10})')
#http://stackoverflow.com/questions/2600790/multiple-levels-of-collection-defaultdict-in-python
#make a nested defaultdict that returns another defaultdict that returns a dict

file_dict = defaultdict(lambda: defaultdict(dict))

root_dir='/newtera/tera/phil/nchaparr/python/Plotting'
varnames=['theta_bar','press','wvelthetapert']
run_name_list = ["Nov302013","Dec142013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]

#
# make a triply nested dict of all vertial profile files
# index is [casename][varname][time]
#
for case in run_name_list:
    for var_name in varnames:
        glob_expr='{}/{}/data/{}*'.format(root_dir,case,var_name)
        files=glob.glob(glob_expr)
        filenames={}
        for file_name in files:
            dirs,name=os.path.split(file_name)
            varname,the_time=expr.match(name).groups(1)
            filenames[int(the_time)]=file_name
        file_dict[case][var_name]=filenames
        print(case,var_name,len(file_dict[case][var_name]))

#all runs have the same heights, choose one at random
height_file='/newtera/tera/phil/nchaparr/python/Plotting/Mar52014/data/heights0000013200'
heights=np.genfromtxt(height_file)
num_levels,=heights.shape

time_full=np.linspace(600,28800,48)
time_nov302013=np.linspace(900,28800,32)

#
# create dataframes for each variable for each case
# first column is height, subsequent columns are times
#
for case in run_name_list:
    for var_name in varnames:
        df=pd.DataFrame(heights,columns=['height'])
        if case == "Nov302013":
            time_list=time_nov302013
        else:
            time_list=time_full
        for the_time in time_list:
            df[int(the_time)]=np.genfromtxt(file_dict[case][var_name][the_time])
        file_dict[case][var_name]['df']=df
#
# now do derived variables for each case and time  rhow and dthetadz
#
varnames.extend(['rhow','dthetadz'])            

for case in run_name_list:
    if case == "Nov302013":
        time_list=time_nov302013
    else:
        time_list=time_full
    df_rhow=pd.DataFrame(heights,columns=['height'])
    df_dthetadz=pd.DataFrame(heights,columns=['height'])
    for the_time in time_list:
        press=file_dict[case]['press']['df'][the_time]
        height=file_dict[case]['press']['df']['height']
        theta=file_dict[case]['theta_bar']['df'][the_time]
        rhow = nc.calc_rhow(press, height, theta[0])
        df_rhow[the_time]=rhow
        dheight = np.diff(height)
        dtheta = np.diff(theta)      
        dthetadz = np.divide(dtheta, dheight)        
        element0 = np.array([0])
        dthetadz=np.hstack((element0, 1.0*dthetadz)) #*1.0/gamma
        df_dthetadz[the_time]=dthetadz
    new_var=file_dict[case]['rhow']
    new_var['df']=df_rhow
    new_var=file_dict[case]['dthetadz']
    new_var['df']=df_dthetadz


#now do the layer time-series variables
layervars=['invrinos','AvProfLims']

thickness_columns=['h0','h','h1','zf0','zf','zf1','deltatheta','mltheta']
rino_columns=['rino', 'invrino', 'wstar', 'S', 'tau', 'mltheta', 'deltatheta', 
              'pi3', 'pi4', 'thetastar']
layer_columns=dict(invrinos=rino_columns,AvProfLims=thickness_columns)
varnames.extend(layervars)            


with pd.HDFStore('new.h5','r') as store:
for case in run_name_list:
    for var_name in layervars:
        colnames=layer_columns[var_name]
        glob_expr='{}/{}/data/{}'.format(root_dir,case,var_name)
        glob_out=glob.glob(glob_expr)
        the_var=np.genfromtxt(glob_out[0])
        if the_var.shape[0]==32:
            big_array=np.empty([48,len(colnames)],dtype=np.float64)
            big_array[:,:]=np.nan
            big_array[0:32,:]=the_var[:,:]
            the_var=big_array[:,:]
        new_var=file_dict[case][var_name]
        new_var['df']=pd.DataFrame.from_records(the_var,columns=colnames)
        node_name='{}/AvProfLims'.format(

#
# write the case information(stability, etc.)
#
with pd.HDFStore('paper_table.h5','r') as store:
     print(store.keys())
     df_cases=store.get('cases')

h5file='all_profiles.h5'
with pd.HDFStore(h5file,'w') as store:
    #loop over text files files
    store.put('df_overview',df_cases,format='table')
    #
    # write the vertical profiles
    # 
    for case in run_name_list:
        for var_name in varnames:
            node_name='/{}/{}'.format(case,var_name)
            store.put(node_name,file_dict[case][var_name]['df'],format='table')            

#
# hdf doesn't know how to store a list of strings, turn it into a single string
# using repr
#
group_attributes={}
history='written 2015/8/5 by plot_theta_profs.py  9542a821e'        
group_attributes['/']=dict(time600=time_full,time900=time_nov302013,history=history,
                           case_list=repr(run_name_list),varnames=repr(varnames))

#
# reopen the file with h5py so we can write attributes to the root node '/'
# (can't figure out how to do this with HDFStore)
#
rootname='/'
with h5py.File(h5file,'a') as f:
    group=f[rootname]
    for key,value in group_attributes[rootname].items():
        print(key,value,type(key),type(value))
        group.attrs[key]=value
