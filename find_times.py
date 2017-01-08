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

h5file='paper_table.h5'
#
# read in the dataframe with
# ['name', 'L0', 'Nperiod', 'N', 'fluxes', 'gammas']
# for each of the 7 cases
#
with pd.HDFStore(h5file,'r') as store:
    print(store.keys())
    nodename='cases'
    df_overview=store.get(nodename)
    Nvals=df_overview['N']
    names=df_overview['name']
    L0=df_overview['L0']
    N_dict={k:v for k,v in zip(names,Nvals)}
    L0_dict={k:v for k,v in zip(names,L0)}
    L0_legend={k:'{:2d}'.format(int(np.round(v,decimals=0))) for k,v in L0_dict.items()}


h5file='good.h5'
df_dict=od()
with pd.HDFStore(h5file,'r') as store:
    df_overview=store.get('/df_overview')
    varnames=list(store.get('/var_names'))
    case_names=list(store.get('/date_names'))
    time600=np.array(store.get('/time600'))
    time900=np.array(store.get('/time900'))
    for case in case_names:
        for name in ['AvProfLims']:
            nodename='{}/{}'.format(case,name)
            df_dict[case]=store.get(nodename)


case_names=df_overview['name']
df_dict={}
for case in case_names:
    if case == 'Nov302013':
        times=time900
    else:
        times=time600
    N=N_dict[case]
    nd_time=times*N
    L0_val=L0_dict[case]
    zenc=find_zenc(times,N,L0_val)
    zenc_L0=zenc/L0_val
    case_dict=dict(zenc=zenc,zenc_L0=zenc_L0,nd_time=nd_time,time_sec=times)
    df_dict[case]=pd.DataFrame(case_dict)
    
target_time=22.5
plt.close('all')
fig,ax = plt.subplots(1,1)
outfile='time_tables.h5'
with pd.HDFStore(outfile,'w') as out:
    for case in case_names:
        df=df_dict[case]
        out.put(case,df)
        zenc_L0=df['zenc_L0'].values
        nd_times=df['nd_time'].values
        index=np.searchsorted(zenc_L0,[target_time])
        print(case,len(zenc_L0),index,nd_times[index])
        ax.plot(zenc_L0,label=case)
ax.legend(loc='best')
plt.show()
        
with h5py.File(outfile,'a') as out:
    out.attrs['history']='written by find_times.py'
print('done')

    
    

    
    
    
    
