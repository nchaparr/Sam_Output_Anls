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
import h5py
from collections import OrderedDict as od
import ast

h5file='paper_table.h5'
with pd.HDFStore(h5file,'r') as store:
    nodename='cases'
    df_overview=store.get(nodename)
    Nvals=df_overview['N']
    names=df_overview['name']
    L0=df_overview['L0']
    N_dict={k:v for k,v in zip(names,Nvals)}
    L0_dict={k:v for k,v in zip(names,L0)}
    L0_legend={k:'{:2d}'.format(int(np.round(v,decimals=0))) for k,v in L0_dict.items()}

time_600=np.linspace(600,28800,48)
time_900=np.linspace(900,28800,32)
varnames=['theta_bar','press','wvelthetapert']
case_list = ["Nov302013","Dec142013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]
case_list.sort()

root_dir="/newtera/tera/phil/nchaparr/python/Plotting"
prof_names=['heights','theta_bar','press','wvelthetapert']
prof_dict=od()
avprof_cols=['h0','h','h1','zf0','zf','zf1','deltatheta','mltheta']
for date in case_list:
    if date == 'Nov302013':
        times=time900
    else:
        times=time600
    dir_path='{}/{}/data'.format(root_dir,date)
    full_path='{:s}/AvProfLims'.format(dir_path)
    numbers=np.genfromtxt(full_path)

    columns=times
    prof_dict[date,'AvProfLims']=numbers
    prof_dict[date,'AvProfLims2']=pd.DataFrame(numbers,columns=avprof_cols)
    sfc_flux=df_overview[df_overview['name']==date]['fluxes']  #W/m^2
    gamma=float(df_overview[df_overview['name']==date]['gammas']/1.e3)  #K/m
    for var in prof_names:
        for the_time in times:
            the_var='{:s}{:010d}'.format(var,int(the_time))
            full_path='{}/{}'.format(dir_path,the_var)
            numbers=np.genfromtxt(full_path)
            prof_dict[date,var,the_time]=numbers
    for index,the_time in enumerate(times):
        theta=prof_dict[date,'theta_bar',the_time]
        press = prof_dict[date,'press',the_time]
        height=prof_dict[date,'heights',the_time]
        rhow = nc.calc_rhow(press, height, theta[0])
        wvelthetapert = prof_dict[date,'wvelthetapert',the_time]
        wvelthetapert[0] = np.nan
        dthetadz = np.diff(theta)/np.diff(height)
        dthetadz=np.hstack((dthetadz, [0]))
        the_h=prof_dict[date,'AvProfLims2'].loc[index]['h']
        scaled_height = height/the_h
        fluxes = wvelthetapert*rhow*1004.0/sfc_flx
        prof_dict[date,'scaled_flux',the_time]=fluxes
        prof_dict[date,'dthetadz',the_time]=dthetadz
        prof_dict[date,'scaled_dtheta',the_time]=dthetadz/gamma
        prof_dict[date,'scaled_height',the_time]=scaled_height
        prof_dict[date,'rhow',the_time]=rhow
        prof_dict[date,'h_val',the_time]=the_h

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

for date in case_list:

    Fig1,axes = plt.subplots(1,3)
    Fig1.suptitle(date)
    for the_ax in axes:
        the_ax.set_ylim(0.1,1.4)

    Ax,Ax1,Ax2=axes

    Ax.set_xlabel(r"$\overline{\theta}$", fontsize=20)
    Ax.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
    Ax.set_xlim(300, 312)

    Ax1.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
    Ax1.set_xticks([.02, 1])
    Ax2.set_xlabel(r"$\frac{\overline{w^{'}\theta^{'}}}{\overline{w^{'}\theta^{'}}_{0}}$", fontsize=20)

    if date == 'Nov302013':
        times=time900
    else:
        times=time600

    for index,the_time in enumerate(times):
        theta=prof_dict[date,'theta_bar',the_time]
        press = prof_dict[date,'press',the_time]
        height=prof_dict[date,'heights',the_time]
        rhow = prof_dict[date,'rhow',the_time]
        wvelthetapert = prof_dict[date,'wvelthetapert',the_time]
        dthetadz = prof_dict[date,'scaled_dtheta',the_time]
        scaled_height = prof_dict[date,'scaled_height',the_time]
        fluxes = prof_dict[date,'scaled_flux',the_time]
        if np.mod(the_time,3600) == 0:
            Ax.plot(theta, scaled_height, '-')
            Ax1.plot(dthetadz,scaled_height, '-')
            Ax2.plot(fluxes, scaled_height, '-')

plt.show()






    
    
