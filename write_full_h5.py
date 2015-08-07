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



def find_zenc(time_sec,N,L0_val):
    zenc=L0_val*np.sqrt(2*time_sec*N)
    return zenc


#
# read the case information(stability, etc.) written by gm_numbers.py
#
with pd.HDFStore('paper_table.h5','r') as store:
     df_overview=store.get('cases')

Nvals=df_overview['N']
names=df_overview['name']
L0=df_overview['L0']
N_dict={k:v for k,v in zip(names,Nvals)}
L0_dict={k:v for k,v in zip(names,L0)}

#split 'theta_bar0000003420' into two parts
expr=re.compile('(.*)?(\d{10,10})')
#http://stackoverflow.com/questions/2600790/multiple-levels-of-collection-defaultdict-in-python
#make a nested defaultdict that returns another defaultdict that returns a dict

file_dict = defaultdict(lambda: defaultdict(dict))

root_dir='/newtera/tera/phil/nchaparr/python/Plotting'
varnames=['theta_bar','press','wvelthetapert']
run_name_list = ["Nov302013","Dec142013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]
run_name_list.sort()

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

#all runs have the same heights, choose one at random
height_file='/newtera/tera/phil/nchaparr/python/Plotting/Mar52014/data/heights0000013200'
heights=np.genfromtxt(height_file)
num_levels,=heights.shape

time_full=np.linspace(600,28800,48)
time_nov302013=np.linspace(900,28800,32)

time_dict={32:time_nov302013,48:time_full}

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
# using the dataframes we've just created
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

for case in run_name_list:
    for var_name in layervars:
        glob_expr='{}/{}/data/{}'.format(root_dir,case,var_name)
        files=glob.glob(glob_expr)
        the_array=np.genfromtxt(files[0])
        df=pd.DataFrame(the_array,columns=layer_columns[var_name])
        if var_name == 'AvProfLims':
            print('all: ',files[0],the_array[0,0])
        df['time_secs']=time_dict[len(df)]
        df['time_nd']=df['time_secs']*N_dict[case]
        df['zenc']=find_zenc(df['time_secs'],N_dict[case],L0_dict[case])
        file_dict[case][var_name]['df']=df

varnames.extend(layervars)            

h5file='all_profiles.h5'
with pd.HDFStore(h5file,'w') as store:
    store.put('df_overview',df_overview,format='table')
    #
    # write  vertical profiles and layer variables
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
        group.attrs[key]=value
