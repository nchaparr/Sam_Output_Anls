import pandas as pd
import h5py
import ast
import numpy as np
from collections import defaultdict,OrderedDict

def find_zenc(time_sec,N,L0_val):
    zenc=L0_val*np.sqrt(2*time_sec*N)
    return zenc

h5file='all_profiles.h5'

#
# get the root attributes
#
with h5py.File(h5file,'r') as f:
    time600=f.attrs['time600']
    time900=f.attrs['time900']
    #
    # turn repr strings into python list objects
    #
    varnames=ast.literal_eval(f.attrs['varnames'])
    case_list=ast.literal_eval(f.attrs['case_list'])

#
# read in the case variables
#

with pd.HDFStore(h5file,'r') as store:
    df_overview=store.get('/df_overview')
#
# get N and L0 for each case
#
Nvals=df_overview['N']
names=df_overview['name']
L0=df_overview['L0']
N_dict={k:v for k,v in zip(names,Nvals)}
L0_dict={k:v for k,v in zip(names,L0)}
L0_legend={k:'{:2d}'.format(int(np.round(v,decimals=0))) for k,v in L0_dict.items()}

#
# for each case, find the time for N*time==200
#
index_dict=OrderedDict()
times_dict=OrderedDict()

#
# find index for each case where non-dimensional time= 250
#
for case in case_list:
    N=N_dict[case]
    if case == 'Nov302013':
        times=time900
        Ntimes=N*time900
    else:
        times=time600
        Ntimes=N*time600
    delta=np.diff(Ntimes)[0]
    index=np.where(np.abs(Ntimes - 250) < delta/2.)[0][0]
    index_dict[case]=index
    times_dict[case]=times[index]
#
# sort cases by L0
#
def case_sort(casename):
    return L0_dict[casename]

case_list.sort(key=case_sort)
#
# now get the profiles  of gradient and flux and the h0,h1 limits
# from dataframes  /case/dthetadz, /case/wvelthetapert /case/AvProfLims
#
h0_dict=OrderedDict()
zf0_dict=OrderedDict()
profile_dict=defaultdict(dict)
h0_fulldict=OrderedDict()
with pd.HDFStore(h5file,'r') as store:
    profile_dict[case]={}
    for case in case_list:
        if case == 'Nov302013':
            times=time900
        else:
            times=time600
        Ntimes=N_dict[case]*times
        nodename='/{}/dthetadz'.format(case)
        df_gradient=store.get(nodename)
        profile_dict[case]['gr']=df_gradient[times_dict[case]]
        nodename='/{}/wvelthetapert'.format(case)
        df_flux=store.get(nodename)
        nodename='/{}/dthetadz'.format(case)
        df_gradient=store.get(nodename)
        profile_dict[case]['flux']=df_flux[times_dict[case]]
        profile_dict[case]['gradient']=df_gradient[times_dict[case]]
        heights=df_flux['height']
        nodename='/{}/AvProfLims'.format(case)
        df_lims=store.get(nodename)
        zenc=find_zenc(times,N_dict[case],L0_dict[case])
        h=df_lims['h'][index_dict[case]]
        h0=df_lims['h0'][index_dict[case]]
        zf0=df_lims['zf0'][index_dict[case]]
        h0_dict[case]=h0/h
        zf0_dict[case]=zf0/h
        h0_fulldict[case]=(Ntimes,df_lims['h0']/zenc)
        profile_dict[case]['nd_heights']=heights/h

from matplotlib import pyplot as plt
plt.close('all')
fig1,ax1=plt.subplots(1,1)
for case in case_list:
    ax1.plot(profile_dict[case]['flux'],profile_dict[case]['nd_heights'],label=L0_legend[case])
    ax1.plot(-0.01,zf0_dict[case],'k*')
ax1.plot([0,0],[1.2,0],'k',linewidth=2)
ax1.set_ylim([0,1.2])
ax1.legend(loc='best')

fig2,ax2=plt.subplots(1,1)
for case in case_list:
    ax2.plot(profile_dict[case]['gradient'],profile_dict[case]['nd_heights'],label=L0_legend[case])
    ax2.plot(-0.0005,h0_dict[case],'k*')
ax2.plot([0,0],[1.2,0],'k',linewidth=2)
ax2.set_ylim([0,1.2])
ax2.set_xlim([-0.001,0.005])
ax2.legend(loc='lower right')


fig3,ax3=plt.subplots(1,1)
for case in case_list:
    ax3.plot(h0_fulldict[case][0],h0_fulldict[case][1],label=L0_legend[case])
ax3.legend(loc='best')
ax3.set_ylim([0.6,1.4])
ax3.set_title('h0 non dimensional')

    
fig4,ax4=plt.subplots(1,1)
with pd.HDFStore(h5file,'r') as store:
    for case in case_list:
        nodename='/{}/AvProfLims'.format(case)
        df_all=store.get(nodename)
        ax4.plot(df_all['time_secs'],df_all['zenc'],label=L0_legend[case])
    ax4.legend(loc='best')
    ax4.set_title('zenc')


plt.show()

