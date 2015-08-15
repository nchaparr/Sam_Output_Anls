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

def get_fields(date,filelist,timelist,varname):
    itervars=zip(timelist,filelist)
    out_dict=od()
    for the_time, filename in itervars:
        out_dict[case,the_time,varname]=np.genfromtxt(filename)
    return out_dict

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
prof_names=['theta_bar','press','wvelthetapert']
for date in case_list:
    if date == 'Nov302013':
        times=time900
    else:
        times=time600
    dir_path='{}/{}/data'.format(root_dir,date)
    for var in prof_names:
        for the_time in times:
            the_var='{:s}{:010d}'.format(var,int(the_time))
            full_path='{}/{}'.format(dir_path,the_var)
            numbers=np.genfromtxt(full_path)

date = "Mar12014"
sfc_flx = 60
gamma = .01

plt.close('all')
Fig1 = plt.figure(1)
Fig1.clf()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


Ax = Fig1.add_subplot(131)
#Ax.set_title( r'$\theta$', fontsize=20)
#Ax.set_title( r'$\frac{\partial \theta}{\partial z}$', fontsize=20)
#Ax.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
Ax.set_xlabel(r"$\overline{\theta}$", fontsize=20)
Ax.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#Ax.set_ylabel(r"$z$", fontsize=20)
plt.xlim(300, 312)
#plt.ylim(100, 1600)
plt.ylim(0.1, 1.4)


Ax1 = Fig1.add_subplot(132)
#Ax1.set_title( r'$Scaled \ \frac{\partial \theta}{\partial z}$', fontsize=20)
#Ax1.set_title( r'$\frac{\partial \theta}{\partial z}$', fontsize=20)
Ax1.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
#Ax1.set_xlabel(r"$\frac{\partial \theta}{\partial z}$", fontsize=20)
#Ax1.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#start, end = -.025, .025
start, end = -1, 2.5
Ax1.set_xticks([.02, 1])
#Ax1.set_ylabel(r"$z$", fontsize=20)
#plt.xlim(-.025, .025)
#plt.xlim(-1, 2.5)
#plt.ylim(100, 1600)
plt.ylim(0.1, 1.4)

Ax2 = Fig1.add_subplot(133)
#Ax2.set_title(r"$\overline{w^{'} \theta^{'}}$", fontsize=20)
#Ax2.set_title(r"$Scaled \ \overline{w^{'} \theta^{'}}$", fontsize=20)
#Ax2.set_xlabel(r"$\overline{w^{'}\theta^{'}}$", fontsize=20)
Ax2.set_xlabel(r"$\frac{\overline{w^{'}\theta^{'}}}{\overline{w^{'}\theta^{'}}_{0}}$", fontsize=20)
#start, end = -.06, .14
#start, end = -.4, 1.2
#Ax2.set_xticks(np.arange(start, end, 1.0*(end-start)/5))

#Ax2.set_ylabel(r"$z$", fontsize=20)
#Ax2.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#plt.ylim(100, 1600)
#plt.xlim(-.06, .14)
#plt.xlim(-.4, 1.2)
plt.ylim(0.1, 1.4)
dump_time_list, Times = Make_Timelists(1, 600, 28800)
 
theta_file_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/theta_bar"+ dump_time for dump_time in dump_time_list]
press_file_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/press"+ dump_time for dump_time in dump_time_list]
flux_file_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/wvelthetapert"+ dump_time for dump_time in dump_time_list]
height_file = "/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/heights0000000600"
AvProfVars = np.genfromtxt("/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/AvProfLims")



#loop over text files files
height = np.genfromtxt(height_file)
df_fluxprof=pd.DataFrame(height,columns=['height'])
df_rhowprof=pd.DataFrame(height,columns=['height'])
h0=[]
time_list=[]
for i in range(len(theta_file_list)):
    print(flux_file_list[i])
    theta = np.genfromtxt(theta_file_list[i])
    press = np.genfromtxt(press_file_list[i])
    rhow = nc.calc_rhow(press, height, theta[0])
    wvelthetapert = np.genfromtxt(flux_file_list[i])
    wvelthetapert[0] = np.nan
    the_time=int(dump_time_list[i])
    time_list.append(the_time)
    col_name='hf_{}'.format(the_time)
    df_fluxprof[the_time]=wvelthetapert
    #Now for the gradients
    dheight = np.diff(height)
    dtheta = np.diff(theta)      
    dthetadz = np.divide(dtheta, dheight)

    element0 = np.array([0])
    dthetadz=np.hstack((dthetadz, element0))

    #only need up to 2500meters
    top_index = np.where(abs(2545 - height) < 40.)[0][0]

    #where gradient is max, and flux is min
    #print AvProfVars[:,1].shape, height.shape
    h0.append(AvProfVars[i,1])
    scaled_height = [1.0*h/AvProfVars[i,1] for h in height]
    df_rhowprof[the_time]=rhow

    fluxes = np.multiply(wvelthetapert, rhow)*1004.0/sfc_flx

    if np.mod(i+1, 6) == 0:
    #if i > 14 and i < 21:

        fluxes[0] = np.nan
        zeros = np.zeros_like(height)

        Ax.plot(theta, scaled_height, '-') #, label = str(Times[i])+'hrs'

        Ax1.plot(1.0*dthetadz/gamma, scaled_height, '-', label = str(Times[i])+'hrs')

        Ax2.plot(fluxes, scaled_height, '-', label = str(Times[i])+'hrs')    

h5file='mar12014.h5'
with pd.HDFStore(h5file,'w') as store:
    nodename='/Mar12014/flux_prof'
    store.put(nodename,df_fluxprof,format='table')
    nodename='/Mar12014/rhow_prof'
    store.put(nodename,df_rhowprof,format='table')
    h0_array=np.array([time_list,h0])
    df_h0=pd.DataFrame(h0_array.T,columns=['times','h0'])
    nodename='/Mar12014/h0'
    store.put(nodename,df_h0,format='table')

rootname='/'
with h5py.File(h5file,'a') as f:
    group=f[rootname]
    group.attrs['comments']='flux units K (m/s), rhow units kg/m^3'
    
array = np.genfromtxt('/newtera/tera/phil/nchaparr/python/Pert_Files/snd')
    
height_0 = array[:, 0]
theta_0 = array[:, 1]
f=interp1d(height_0, theta_0)

#Now plot inital sounding
top_index = np.where(height <= 2500)[0][-1]
theta_0 = f(height[0:top_index])
dtheta0 = np.diff(theta_0)
dthetadz0 = np.divide(dtheta0, dheight[0:top_index-1])
element0 = np.array([.005])
dthetadz0=np.hstack((element0, dthetadz0))

#Ax1.plot(dthetadz0, scaled_height[0:top_index], '--', label = 'Initial Sounding')
Ax1.plot(zeros, height)#zeros line for reference
#Ax1.plot(gamma, scaled_height)#zeros line for reference
Ax1.plot(zeros+1, height, 'k-')#zeros line for reference
Ax1.plot(zeros+.02, height, 'k-')#zeros line for reference
#Ax2.plot(zeros, height)#zeros line for reference
Ax2.plot(zeros, height)#zeros line for reference 
plt.legend(loc = 'lower right', prop={'size':8})
#Ax2.plot(theta_0, scaled_xheight[0:top_index], '--', label = 'Initial Sounding')#"
#plt.xlim(300, 310)
plt.legend(loc = 'upper right', prop={'size':8})
plt.show()
#Fig1.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/theta_flux_profs.png")





    
    
