import pandas as pd
import h5py
import ast
import numpy as np
from collections import OrderedDict as od


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
profile_dict=od()
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
        nodename='/{}/rhow'.format(case)
        df_rhow=store.get(nodename)
        profile_dict[case,'rhow']=df_rhow
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
# get all the variables we want to plot
# note that every case was run with the same vertical grid
#
heights=profile_dict[case_list[0],'gradient']['height']
for case in compare_cases:
    zenc=profile_dict[case,'zenc']
    for nd_time in target_times:
        index=index_dict[nd_time,case]
        plot_dict[case,nd_time,'df_prof']=pd.DataFrame(heights,columns=['height'])
        plot_dict[case,nd_time,'h']=profile_dict[case,'df_lims']['h'][index]
        plot_dict[case,nd_time,'h0']=profile_dict[case,'df_lims']['h0'][index]
        plot_dict[case,nd_time,'zenc']=zenc[index]
        plot_dict[case,nd_time,'zf0']=profile_dict[case,'df_lims']['zf0'][index]
        the_time=times_dict[nd_time,case]
        plot_dict[case,nd_time,'df_prof']['flux_prof']=profile_dict[case,'flux'][the_time]
        plot_dict[case,nd_time,'df_prof']['grad_prof']=profile_dict[case,'gradient'][the_time]
        plot_dict[case,nd_time,'df_prof']['rhow']=profile_dict[case,'rhow'][the_time]

from matplotlib import pyplot as plt
plt.close('all')

cpd=1004.
for the_time in target_times:
    flux_fig,flux_ax=plt.subplots(1,1)
    grad_fig,grad_ax=plt.subplots(1,1)
    for case in compare_cases:
        sfc_flux=float(df_overview[df_overview['name']==case]['fluxes'])
        df_prof=plot_dict[case,the_time,'df_prof']
        rhow=plot_dict[case,the_time,'df_prof']['rhow']
        df_prof['ndh_flux']=df_prof['height']/plot_dict[case,nd_time,'zf0']
        df_prof['nd_flux']=plot_dict[case,nd_time,'df_prof']['flux_prof']*rhow*cpd/sfc_flux
        df_prof['ndh_h']=df_prof['height']/plot_dict[case,nd_time,'h']
        df_prof.plot(x='nd_flux',y='ndh_flux',label=L0_legend[case],ax=flux_ax,legend=False)
        df_prof.plot(x='grad_prof',y='ndh_h',label=L0_legend[case],ax=grad_ax,legend=False)
    flux_ax.set_title('flux profile at nd time of {}'.format(the_time))
    grad_ax.set_title('gradient profile at nd time of {}'.format(the_time))
    flux_ax.set_ylim(0,1.3)
    grad_ax.set_ylim(0,1.3)
    flux_ax.set_xlim([-0.3,1.])
    grad_ax.set_xlim([-0.001,0.005])
    flux_ax.legend(loc='best')
    grad_ax.legend(loc='best')
    filename='flux_{}.png'.format(the_time)
    flux_fig.savefig(filename)
    filename='grad_{}.png'.format(the_time)
    grad_fig.savefig(filename)
plt.show()



