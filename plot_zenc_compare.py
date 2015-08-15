import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from collections import OrderedDict as od
h5new='good.h5'

#
# get the root attributes
#
df_dict=od()
with pd.HDFStore(h5new,'r') as store:
    df_overview=store.get('/df_overview')
    varnames=list(store.get('/var_names'))
    case_list=list(store.get('/date_names'))
    time600=np.array(store.get('/time600'))
    time900=np.array(store.get('/time900'))
    for case in case_list:
        for name in ['AvProfLims']:
            nodename='{}/{}'.format(case,name)
            df_dict[case,name]=store.get(nodename)

Nvals=df_overview['N']
names=df_overview['name']
L0=df_overview['L0']
N_dict={k:v for k,v in zip(names,Nvals)}
L0_dict={k:v for k,v in zip(names,L0)}
L0_legend={k:'{:2d}'.format(int(np.round(v,decimals=0))) for k,v in L0_dict.items()}

n_colors=len(case_list)
colors=plt.cm.rainbow(np.linspace(0,1,n_colors))

def case_sort(casename):
    return L0_dict[casename]
case_list.sort(key=case_sort)
    
plt.close('all')

cp=1004.

fig=plt.figure()
ax=fig.add_axes([0.11,0.1,0.73,0.8])

for i,case in enumerate(case_list):
    df_lims=df_dict[case,'AvProfLims']
    label='{}_h0'.format(L0_legend[case])
    df_lims.plot(x='Ntimes',y='h0',linestyle='None',marker='s',color=colors[i],ax=ax,label=label,legend=False)
    label='{}_zenc'.format(L0_legend[case])
    df_lims.plot(x='Ntimes',y='zenc',linestyle='-',color=colors[i],ax=ax,label=label,legend=False)
ax.legend(numpoints=1,bbox_to_anchor=(1.,0.9),bbox_transform=fig.transFigure,prop={'size':11})
ax.set_ylim([0,1800])
ax.set_ylabel('height (m)')
ax.set_title('h0 vs. zenc for all profiles')
fig.savefig('h0_time.png')

    
fig=plt.figure()
ax=fig.add_axes([0.11,0.1,0.73,0.8])

for i,case in enumerate(case_list):
    df_lims=df_dict[case,'AvProfLims']
    label='{}_h'.format(L0_legend[case])
    df_lims.plot(x='Ntimes',y='h',linestyle='None',marker='s',color=colors[i],ax=ax,label=label,legend=False)
    label='{}_zenc'.format(L0_legend[case])
    df_lims.plot(x='Ntimes',y='zenc',linestyle='-',color=colors[i],ax=ax,label=label,legend=False)
ax.legend(numpoints=1,bbox_to_anchor=(1.,0.9),bbox_transform=fig.transFigure,prop={'size':11})
ax.set_ylim([0,1800])
ax.set_ylabel('height (m)')
ax.set_title('h vs. zenc for all profiles')
fig.savefig('h_time.png')


fig=plt.figure()
ax=fig.add_axes([0.11,0.1,0.73,0.8])

for i,case in enumerate(case_list):
    df_lims=df_dict[case,'AvProfLims']
    label='{}_h1'.format(L0_legend[case])
    df_lims.plot(x='Ntimes',y='h1',linestyle='None',marker='s',color=colors[i],ax=ax,label=label,legend=False)
    label='{}_zenc'.format(L0_legend[case])
    df_lims.plot(x='Ntimes',y='zenc',linestyle='-',color=colors[i],ax=ax,label=label,legend=False)
ax.legend(numpoints=1,bbox_to_anchor=(1.,0.9),bbox_transform=fig.transFigure,prop={'size':11})
ax.set_ylim([0,1800])
ax.set_ylabel('height (m)')
ax.set_title('h1 vs. zenc for all profiles')
fig.savefig('h1_time.png')

    
    
plt.show()

