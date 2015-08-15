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

"""

   For plotting the (scaled) temperature gradient and flux profiles.

"""

h5file='all_profiles.h5'
with pd.HDFStore(h5file,'r') as store:
    nodename='/df_overview'
    df_overview=store.get(nodename)
    Nvals=df_overview['N']
    names=df_overview['name']
    L0=df_overview['L0']
    N_dict={k:v for k,v in zip(names,Nvals)}
    L0_dict={k:v for k,v in zip(names,L0)}
    L0_legend={k:'{:2d}'.format(int(np.round(v,decimals=0))) for k,v in L0_dict.items()}

with h5py.File(h5file,'r') as f:
    time600=f.attrs['time600']
    time900=f.attrs['time900']
    #
    # turn repr strings into python list objects
    #
    varnames=ast.literal_eval(f.attrs['varnames'])
    case_list=ast.literal_eval(f.attrs['case_list'])

root_dir="/newtera/tera/phil/nchaparr/python/Plotting"
prof_names=['heights','theta_bar','press','wvelthetapert']
prof_dict=od()
for date in case_list:
    if date == 'Nov302013':
        times=time900
    else:
        times=time600
    dir_path='{}/{}/data'.format(root_dir,date)
    full_path='{:s}/AvProfLims'.format(dir_path)
    numbers=np.genfromtxt(full_path)
    prof_dict[date,'AvProfLims']=numbers
    for the_time in times:
        for var in prof_names:
            the_var='{:s}{:010d}'.format(var,int(the_time))
            full_path='{}/{}'.format(dir_path,the_var)
            numbers=np.genfromtxt(full_path)
            prof_dict[date,var,the_time]=numbers

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Fig1,axes = plt.subplots(1,3)
for the_ax in axes:
    the_ax.set_ylim(0.1,1.4)

Ax,Ax1,Ax2=axes

Ax.set_xlabel(r"$\overline{\theta}$", fontsize=20)
Ax.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
Ax.set_xlim(300, 312)

Ax1.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
start, end = -1, 2.5
Ax1.set_xticks([.02, 1])
Ax2.set_xlabel(r"$\frac{\overline{w^{'}\theta^{'}}}{\overline{w^{'}\theta^{'}}_{0}}$", fontsize=20)

date = "Mar12014"
sfc_flux=df_overview[df_overview['name']==date]['fluxes']  #W/m^2
gamma=float(df_overview[df_overview['name']==date]['gammas']/1.e3)  #K/m


for index,the_time in enumerate(time600):
    theta=prof_dict[date,'theta_bar',the_time]
    press = prof_dict[date,'press',the_time]
    height=prof_dict[date,'heights',the_time]
    rhow = nc.calc_rhow(press, height, theta[0])
    wvelthetapert = prof_dict[date,'wvelthetapert',the_time]
    wvelthetapert[0] = np.nan
    #Now for the gradients
    dheight = np.diff(height)
    dtheta = np.diff(theta)      
    dthetadz = dtheta/dheight

    element0 = np.array([0])
    dthetadz=np.hstack((dthetadz, element0))

    the_h=prof_dict[date,'AvProfLims'][index,1]
    scaled_height = [1.0*h/the_h for h in height]
    fluxes = wvelthetapert*rhow*1004.0/sfc_flx
    if np.mod(the_time,3600) == 0:
        print(the_time,the_h)
        Ax.plot(theta, scaled_height, '-') #, label = str(Times[i])+'hrs'
        Ax1.plot(1.0*dthetadz/gamma, scaled_height, '-', label = str(Times[i])+'hrs')
        Ax2.plot(fluxes, scaled_height, '-', label = str(Times[i])+'hrs')    

plt.show()






    
    
