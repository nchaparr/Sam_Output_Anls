import pandas as pd
import h5py
import ast
import numpy as np
from collections import OrderedDict as od
from ordered_defaultdict import OrderedDefaultdict as od_def
import pdb


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
index_dict=od()
times_dict=od()
ntimes_dict=od()
#
# find index for each case where non-dimensional time= 250 and 500
#
target_times=[250,500]
for nd_time in target_times:
    for case in case_list:
        N=N_dict[case]
        if case == 'Nov302013':
            times=time900
            Ntimes=N*time900
        else:
            times=time600
            Ntimes=N*time600
            delta=np.diff(Ntimes)[0]
            try:
                index=np.where(np.abs(Ntimes - nd_time) < delta/2.)[0][0]
            except IndexError:
                continue
            index_dict[(nd_time,case)]=index
            times_dict[(nd_time,case)]=times[index]
            ntimes_dict[(nd_time,case)]=Ntimes[index]
#
# these are the cases that have all target times
#
compare_cases=[item[1] for item in times_dict.keys() if item[0]==500]
l0_cases=[L0_dict[case] for case in compare_cases]
#
# sort cases by L0
#
def case_sort(casename):
    return L0_dict[casename]

compare_cases.sort(key=case_sort)
#
# now get the profiles  of gradient and flux and the h0,h1 limits
# from dataframes  /case/dthetadz, /case/wvelthetapert /case/AvProfLims
#
h0_point_dict=od()
zf0_point_dict=od()
profile_dict=od_def(od)
h0_fulldict=od()
#
# read in all profiles for all cases
#
with pd.HDFStore(h5file,'r') as store:
    for case in case_list:
        if case == 'Nov302013':
            times=time900
        else:
            times=time600
        Ntimes=N_dict[case]*times
        nodename='/{}/dthetadz'.format(case)
        df_gradient=store.get(nodename)
        profile_dict[case,'gradient']=df_gradient
        nodename='/{}/wvelthetapert'.format(case)
        df_flux=store.get(nodename)
        profile_dict[case,'flux']=df_flux
        nodename='/{}/dthetadz'.format(case)
        df_gradient=store.get(nodename)
        nodename='/{}/AvProfLims'.format(case)
        df_lims=store.get(nodename)
        profile_dict[case,'df_lims']=df_lims
        print('input: ',case,len(df_lims))
        zenc=find_zenc(times,N_dict[case],L0_dict[case])
        profile_dict[case,'zenc']=zenc
#
# extract vertical profiles of flux and gradient at the
# target times for plotting
#
plot_dict=od()
#
# only one set of dimensional heights
#
heights=profile_dict[case_list[0],'gradient']['height']
for case in compare_cases:
    zenc=profile_dict[case,'zenc']
    for nd_time in target_times:
        index=index_dict[nd_time,case]
        print(nd_time,case,index)
        plot_dict[case,nd_time,'nd_heights']=heights/zenc[index]
        key=(nd_time,case)
        h=profile_dict[case,'df_lims']['h'][index]
        plot_dict[case,nd_time,'h_nd']=h/zenc[index]
        h0=profile_dict[case,'df_lims']['h0'][index]
        zf0=profile_dict[case,'df_lims']['zf0'][index]
        plot_dict[case,nd_time,'h0_nd']=h0/zenc[index]
        plot_dict[case,nd_time,'zf0_nd']=zf0/zenc[index]
        the_time=times_dict[nd_time,case]
        plot_dict[case,nd_time,'flux_prof']=profile_dict[case,'flux'][the_time]
        plot_dict[case,nd_time,'grad_prof']=profile_dict[case,'gradient'][the_time]

from matplotlib import pyplot as plt
plt.close('all')

for the_time in target_times:
    fig1,ax1=plt.subplots(1,1)
    fig2,ax2=plt.subplots(1,1)
    for case in compare_cases:
        x1=plot_dict[case,the_time,'flux_prof']
        x2=plot_dict[case,the_time,'grad_prof']
        y=plot_dict[case,the_time,'nd_heights']
        ax1.plot(x1,y,label=L0_legend[case])
        ax2.plot(x2,y,label=L0_legend[case])
        x=0.
        y1=plot_dict[case,the_time,'zf0_nd']
        y2=plot_dict[case,the_time,'h0_nd']
        ax1.plot(x,y1,'k*')
        ax2.plot(x,y2,'k*')
    ax1.set_ylim(0,1.3)
    ax2.set_ylim(0,1.3)
    ax2.set_xlim([-0.001,0.005])
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax1.set_title('flux profile at nd time of {}'.format(the_time))
    ax2.set_title('gradient profile at nd time of {}'.format(the_time))
plt.show()

#     ax1.plot(-0.01,zf0_point_dict[case],'k*')
# ax1.plot([0,0],[1.2,0],'k',linewidth=2)
# ax1.set_ylim([0,1.2])
# ax1.legend(loc='best')


# fig2,ax2=plt.subplots(1,1)
# for case in case_list:
#     ax2.plot(profile_dict[case]['gradient'],profile_dict[case]['nd_heights'],label=L0_legend[case])
#     ax2.plot(-0.0005,h0_point_dict[case],'k*')
# ax2.plot([0,0],[1.2,0],'k',linewidth=2)
# ax2.set_ylim([0,1.2])
# ax2.set_xlim([-0.001,0.005])
# ax2.legend(loc='lower right')


# fig3,ax3=plt.subplots(1,1)
# for case in case_list:
#     ax3.plot(h0_fulldict[case][0],h0_fulldict[case][1],label=L0_legend[case])
# ax3.legend(loc='best')
# ax3.set_ylim([0.6,1.4])

# ax3.set_title('h0 non dimensional')

    
# fig4,ax4=plt.subplots(1,1)
# with pd.HDFStore(h5file,'r') as store:
#     for case in case_list:
#         nodename='/{}/AvProfLims'.format(case)
#         df_all=store.get(nodename)
#         ax4.plot(df_all['time_secs'],df_all['zenc'],label=L0_legend[case])
#     ax4.legend(loc='best')
#     ax4.set_title('zenc')


# plt.show()

