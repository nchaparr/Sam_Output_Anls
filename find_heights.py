import h5py
from enum import IntEnum
from collections import OrderedDict as od
import numpy as np
from matplotlib import pyplot as plt

#https://github.com/nchaparr/Sam_Output_Anls/blob/master/Readme.txt#L45    
class Heights(IntEnum):
    zg0 = 0
    zg  = 1
    zg1 = 2
    zf0 = 3
    zf =  4
    zf1 = 5
    deltheta = 6
    mltheta=7
    invrinos=8

#
# get all 6 heigths
#
columns=list(Heights.__members__.keys())[:6]

h5_profiles='profiles_flux_quads_260.h5'
with h5py.File(h5_profiles,'r') as infile:
    case_list=list(infile.keys())[:1]
    case=case_list[0]
    height = infile[case]['height'][...]
    hvals = infile[case]['hvals'][...]
    print('case: ',case)
    print('time index: ',infile[case]['time_list'].attrs['time_index'])
    time_index = int(infile[case]['time_list'].attrs['time_index'])
    height_dict=od()
    for col in columns:
        height_dict[col]=hvals[time_index,int(Heights[col])]

index_dict=od()
for key,value in height_dict.items():
    index_dict[key]=np.searchsorted(height,value,side='left')

prof_file='full_out_260.h5'
mean_flux=[]
theta_rms=[]
top_lev=index_dict['zg1']
bot_lev = index_dict['zg0']
levs= list(range(bot_lev-1,top_lev+1))
plot_heights=height[levs]
with h5py.File(prof_file,'r') as flux_in:
    for zlev in levs:
        print('zlev: ',zlev)
        wpert=flux_in[case][str(zlev)]['wpert'][...]
        thetapert=flux_in[case][str(zlev)]['thetapert'][...]
        hit = np.logical_and(thetapert >0, wpert >0)
        theta_uw=thetapert[hit]
        theta_out = np.sqrt(np.mean(theta_uw**2.))
        flux=(wpert*thetapert).mean()
        mean_flux.append(flux)
        theta_rms.append(theta_out)
        
plt.close()

def plot2(prof1,prof2,zlevs,height_dict):
    fig,ax1=plt.subplots(1,1)
    ax2 = ax1.twiny()
    ax1.plot(prof1,zlevs,'b-')
    ax2.plot(prof2,zlevs,'r-')
    for t1 in ax1.get_xticklabels():
        t1.set_color('b')
    for t2 in ax2.get_xticklabels():
        t2.set_color('r')
    for key,hlev in height_dict.items():
        ax1.text(0.,hlev,key)
    return ax1,ax2

zlevs=height[levs]
ax1,ax2 = plot2(mean_flux,theta_rms,zlevs,height_dict)
ax1.set_xlabel('net flux (K m/s)')
ax1.xaxis.label.set_color('b')
ax2.set_xlabel('thetarms (K)')
ax2.xaxis.label.set_color('r')



        
plt.show()
