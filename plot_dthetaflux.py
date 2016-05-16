import numpy as np
import matplotlib.pyplot as plt

import nchap_fun as nc
from matplotlib import rcParams
rcParams.update({'font.size': 10})
import pandas as pd
import h5py
from collections import OrderedDict as od
import plot_gm_numbers,importlib
importlib.reload(plot_gm_numbers)
from plot_gm_numbers import gm_vars
from scipy import signal

def find_zenc(time_sec,N,L0_val):
    zenc=L0_val*np.sqrt(2*time_sec*N)
    return zenc

def moving_average(interval,height, window_size):
    window= np.ones(int(window_size))/float(window_size)
    print('window: ',window)
    return np.convolve(interval, window, 'same')

h5file='paper_table.h5'
#
# read in the dataframe with
# ['name', 'L0', 'Nperiod', 'N', 'fluxes', 'gammas']
# for each of the 7 cases
#
with pd.HDFStore(h5file,'r') as store:
    nodename='cases'
    df_overview=store.get(nodename)
    Nvals=df_overview['N']
    names=df_overview['name']
    L0=df_overview['L0']
    N_dict={k:v for k,v in zip(names,Nvals)}
    L0_dict={k:v for k,v in zip(names,L0)}
    L0_legend={k:'{:2d}'.format(int(np.round(v,decimals=0))) for k,v in L0_dict.items()}

time600=np.linspace(600,28800,48)
time900=np.linspace(900,28800,32)
varnames=['theta_bar','press','wvelthetapert']
case_list = ["Nov302013","Dec142013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]
case_list.sort()

root_dir="/newtera/tera/phil/nchaparr/python/Plotting"
root_dir='/tera/users/nchaparr'
prof_names=['heights','theta_bar','press','wvelthetapert','rootmeanwsq']
prof_dict=od()
#
# read in the horizontal limit files (AvProfLims, time) and
# and the vertical profiles 'heights','theta_bar','press','wvelthetapert'
# (height,time)  and calculate derivied profiles, saving
# the arrays as dataframes
# 
#
avprof_cols=['time', 'h0','h','h1','zf0','zf','zf1','deltatheta','mltheta']
for date in case_list:
    if date == 'Nov302013':
        times=time900
    else:
        times=time600
    dir_path='{}/{}/data'.format(root_dir,date)
    full_path='{:s}/AvProfLims'.format(dir_path)
    numbers=np.genfromtxt(full_path)
    #
    # save the heights to a dataframe
    #
    prof_dict[date,'AvProfLims']=pd.DataFrame(numbers,columns=avprof_cols)
    #
    # add the encroachment height
    #
    zenc=find_zenc(times,N_dict[date],L0_dict[date])
    prof_dict[date,'AvProfLims']['zenc']=zenc
    prof_dict[date,'AvProfLims']['times']=times
    prof_dict[date,'AvProfLims']['Ntimes']=times*N_dict[date]
    sfc_flux=float(df_overview[df_overview['name']==date]['fluxes'])  #W/m^2
    print('here is sfc_flux',sfc_flux)
    gamma=float(df_overview[df_overview['name']==date]['gammas']/1.e3)  #K/m
    the_var='{:s}{:010d}'.format('heights',int(times[0]))
    full_path='{}/{}'.format(dir_path,the_var)
    height=np.genfromtxt(full_path)
    for var in prof_names:
        #
        # initialize the first column sor each new dataframe
        #
        prof_dict[date,var]=pd.DataFrame(height,columns=['height'])
        for the_time in times:
            the_var='{:s}{:010d}'.format(var,int(the_time))
            full_path='{}/{}'.format(dir_path,the_var)
            numbers=np.genfromtxt(full_path)
            prof_dict[date,var][the_time]=numbers
    derived_names=['scaled_flux','dthetadz','scaled_dtheta',
                   'scaled_height','rhow','scaled_wrms','flux_grad','scaled_flux_grad']
    for var in derived_names:
        prof_dict[date,var]=pd.DataFrame(height,columns=['height'])
    cpd=1004.
    for index,the_time in enumerate(times):
        theta=prof_dict[date,'theta_bar'][the_time]
        press = prof_dict[date,'press'][the_time]
        height=prof_dict[date,'heights'][the_time]
        rhow = nc.calc_rhow(press, height, theta[0])
        wvelthetapert = prof_dict[date,'wvelthetapert'][the_time]
        wrms=prof_dict[date,'rootmeanwsq'][the_time]
        wvelthetapert[0] = np.nan
        dthetadz = np.diff(theta)/np.diff(height)
        dthetadz=np.hstack((dthetadz, [0]))
        #flux_smooth=signal.medfilt(wvelthetapert,3)
        flux_smooth=wvelthetapert
        dfluxdz = np.diff(flux_smooth)/np.diff(height)
        dfluxdz=np.hstack((dfluxdz, [0]))
        the_h=prof_dict[date,'AvProfLims'].loc[index]['h']
        scaled_height = height/the_h
        fluxes = wvelthetapert*rhow*cpd/sfc_flux
        theta_0=prof_dict[date,'theta_bar'][the_time][0]
        L0,N,B0,wstar,wstar_gm=gm_vars(sfc_flux,gamma,the_h,theta_0)
        print('L0,wstar',L0,wstar)
        prof_dict[date,'scaled_flux'][the_time]=fluxes
        prof_dict[date,'dthetadz'][the_time]=dthetadz
        prof_dict[date,'scaled_dtheta'][the_time]=dthetadz/gamma
        prof_dict[date,'scaled_height'][the_time]=scaled_height
        prof_dict[date,'scaled_wrms'][the_time]=wrms/wstar
        prof_dict[date,'flux_grad'][the_time]=dfluxdz
        prof_dict[date,'scaled_flux_grad'][the_time]=dfluxdz/gamma
        prof_dict[date,'rhow'][the_time]=rhow
#
# write out all frames to and hdf file using pytables
#
var_list=[]
date_list=[]
all_frames=prof_dict.keys()
var_names=set()
date_names=set()
h5file='good.h5'
with pd.HDFStore(h5file,'w') as store:
    node_name='/df_overview'
    print(df_overview)
    store.put(node_name,df_overview,format='table')
    for date,var in all_frames:
        date_names.add(date)
        var_names.add(var)
        node_name='/{}/{}'.format(date,var)
        store.put(node_name,prof_dict[date,var],format='table')
    var_names=list(var_names)    
    var_names.sort()
    date_names=list(date_names)
    date_names.sort()
    var_ser=pd.Series(var_names)
    date_ser=pd.Series(date_names)
    store.put('/var_names',var_ser)
    store.put('/date_names',date_ser)
    store.put('/time600',pd.Series(time600))
    store.put('/time900',pd.Series(time900))

group_attributes={}
history='written 2015/8/15 by plot_dthetaflux.py  e6dec4235a2c'        
group_attributes['/']=dict(history=history)

rootname='/'
with h5py.File(h5file,'a') as f:
    group=f[rootname]
    for key,value in group_attributes[rootname].items():
        group.attrs[key]=value

#
# now do the plots
#
        

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

for date in case_list:

    Fig1,axes = plt.subplots(1,4)
    title='{} L0={}'.format(date,L0_legend[date])
    Fig1.suptitle(title)
    for the_ax in axes:
        the_ax.set_ylim(0.1,1.4)

    Ax,Ax1,Ax2,Ax3=axes

    Ax.set_xlabel(r"$\overline{\theta}$", fontsize=20)
    Ax.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
    Ax.set_xlim(300, 312)

    Ax1.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
    Ax1.set_xticks([.02, 1])
    Ax2.set_xlabel(r"$\frac{\overline{w^{'}\theta^{'}}}{\overline{w^{'}\theta^{'}}_{0}}$", fontsize=20)
    Ax3.set_xlabel(r"$\frac{\overline{dw^{'}\theta^{'}}}{dz}/\gamma$", fontsize=20)

    if date == 'Nov302013':
        times=time900
    else:
        times=time600

    hour=1
    for index,the_time in enumerate(times):
        scaled_height = prof_dict[date,'scaled_height'][the_time]
        theta=prof_dict[date,'theta_bar'][the_time]
        dthetadz = prof_dict[date,'scaled_dtheta'][the_time]
        fluxes = prof_dict[date,'scaled_flux'][the_time]
        scaled_dflxdz = prof_dict[date,'scaled_flux_grad'][the_time]
        dflxdz=prof_dict[date,'flux_grad'][the_time]
        label='hour {}'.format(hour)
        if np.mod(the_time,3600) == 0:
            Ax.plot(theta, scaled_height, '-',label=label)
            Ax1.plot(dthetadz,scaled_height, '-',label=label)
            Ax2.plot(fluxes, scaled_height, '-',label=label)
            Ax3.plot(scaled_dflxdz,scaled_height,'-',label=label)
            hour+=1
        Ax2.legend(loc='best')
    the_file='L0_{}.png'.format(L0_legend[date])
    Fig1.savefig(the_file)

plt.show()






    
    
