import pandas as pd
import h5py
import ast
import numpy as np
from collections import OrderedDict as od


def find_zenc(time_sec,N,L0_val):
    zenc=L0_val*np.sqrt(2*time_sec*N)
    return zenc

h5file_new='good.h5'
#
# get the root attributes
#
df_dict=od()
with pd.HDFStore(h5file_new,'r') as store:
    df_overview=store.get('/df_overview')
    varnames=list(store.get('/var_names'))
    case_list=list(store.get('/date_names'))
    time600=np.array(store.get('/time600'))
    time900=np.array(store.get('/time900'))
    for case in case_list:
        for name in varnames:
            nodename='{}/{}'.format(case,name)
            df_dict[case,name]=store.get(nodename)
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


from matplotlib import pyplot as plt
cpd=1004.
plt.close('all')

for the_time in target_times:
    flux_fig,flux_ax=plt.subplots(1,1)
    grad_fig,grad_ax=plt.subplots(1,1)
    for case in compare_cases:
        real_time=times_dict[nd_time,case]
        x=df_dict[case,'scaled_flux'][real_time]
        y=df_dict[case,'scaled_height'][real_time]
        flux_ax.plot(x,y,label=L0_legend[case])
        x=df_dict[case,'scaled_dtheta'][real_time]
        grad_ax.plot(x,y,label=L0_legend[case])
    flux_ax.set_title('flux profile at nd time of {}'.format(the_time))
    grad_ax.set_title('gradient profile at nd time of {}'.format(the_time))
    [ax.set_ylim(0,1.3) for ax in [flux_ax,grad_ax]]
    grad_ax.set_ylim(0,1.3)
    flux_ax.set_xlim([-0.3,1.])
    grad_ax.set_xlim([-0.3,2.0])
    flux_ax.legend(loc='best')
    grad_ax.legend(loc='best')
    filename='flux_{}.png'.format(the_time)
    flux_fig.savefig(filename)
    filename='grad_{}.png'.format(the_time)
    grad_fig.savefig(filename)
plt.show()



