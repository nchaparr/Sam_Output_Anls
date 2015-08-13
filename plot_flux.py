import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import h5py,ast

with pd.HDFStore('mar12014.h5','r') as store:
    nodename='/Mar12014/flux_prof'
    df_flux=store.get(nodename)
    nodename='/Mar12014/rhow_prof'
    df_rhow=store.get(nodename)
    nodename='/Mar12014/h0'
    df_h0=store.get(nodename)

with pd.HDFStore('all_profiles.h5','r') as store:
    nodename='/Mar12014//AvProfLims'
    df_mar1_lims=store.get(nodename)
#
# get the root attributes
#
h5file='all_profiles.h5'
with h5py.File(h5file,'r') as f:
    #
    # turn repr strings into python list objects
    #
    varnames=ast.literal_eval(f.attrs['varnames'])
    case_list=ast.literal_eval(f.attrs['case_list'])

with pd.HDFStore(h5file,'r') as store:
    nodename='/df_overview'
    df_overview=store.get(nodename)
    Nvals=df_overview['N']
    names=df_overview['name']
    L0=df_overview['L0']
    N_dict={k:v for k,v in zip(names,Nvals)}
    L0_dict={k:v for k,v in zip(names,L0)}
    L0_legend={k:'{:2d}'.format(int(np.round(v,decimals=0))) for k,v in L0_dict.items()}

hit=df_overview['name']=='Mar12014'
case_values=df_overview[hit]

n_colors=len(case_list)
colors=plt.cm.rainbow(np.linspace(0,1,n_colors))

def case_sort(casename):
    return L0_dict[casename]
case_list.sort(key=case_sort)
    
plt.close('all')

cp=1004.
hit=np.where(df_h0['times']==27600)
h0=float(df_h0.loc[hit]['h0'])
df_flux['height_nd']=df_flux['height']/h0
df_flux['flux_nd']=df_flux[27600]*df_rhow[27600]*cp/float(case_values['fluxes'])

fig,ax=plt.subplots(1,1)
df_flux.plot(y='height_nd',x='flux_nd',ax=ax,legend=False)
ax.set_xlabel('non-dimensional flux')
ax.set_ylabel('height/h0')
ax.set_title('Mar12014 27600 seconds')
ax.set_ylim(0,1.3)
ax.set_xlim(-0.2,1)
fig.savefig('niamh_flux.png')


fig=plt.figure()
ax=fig.add_axes([0.11,0.1,0.73,0.8])

with pd.HDFStore(h5file,'r') as store:
    for i,case in enumerate(case_list):
        nodename='/{}/AvProfLims'.format(case)
        df_lims=store.get(nodename)
        label='{}_h0'.format(L0_legend[case])
        df_lims.plot(x='time_nd',y='h0',linestyle='None',marker='s',color=colors[i],ax=ax,label=label,legend=False)
        label='{}_zenc'.format(L0_legend[case])
        df_lims.plot(x='time_nd',y='zenc',linestyle='-',color=colors[i],ax=ax,label=label,legend=False)
    ax.legend(numpoints=1,bbox_to_anchor=(1.,0.9),bbox_transform=fig.transFigure,prop={'size':11})
    ax.set_ylim([0,1800])
    ax.set_ylabel('height (m)')
    ax.set_title('h0 vs. zenc for all profiles')
fig.savefig('h0_time.png')

    
fig=plt.figure()
ax=fig.add_axes([0.11,0.1,0.73,0.8])

with pd.HDFStore(h5file,'r') as store:
    for i,case in enumerate(case_list):
        nodename='/{}/AvProfLims'.format(case)
        df_lims=store.get(nodename)
        label='{}_h'.format(L0_legend[case])
        df_lims.plot(x='time_nd',y='h',linestyle='None',marker='s',color=colors[i],ax=ax,label=label,legend=False)
        label='{}_zenc'.format(L0_legend[case])
        df_lims.plot(x='time_nd',y='zenc',linestyle='-',color=colors[i],ax=ax,label=label,legend=False)
    ax.legend(numpoints=1,bbox_to_anchor=(1.,0.9),bbox_transform=fig.transFigure,prop={'size':11})
    ax.set_ylim([0,1800])
    ax.set_ylabel('height (m)')
    ax.set_title('h vs. zenc for all profiles')
fig.savefig('h_time.png')


fig=plt.figure()
ax=fig.add_axes([0.11,0.1,0.73,0.8])

with pd.HDFStore(h5file,'r') as store:
    for i,case in enumerate(case_list):
        nodename='/{}/AvProfLims'.format(case)
        df_lims=store.get(nodename)
        label='{}_h1'.format(L0_legend[case])
        df_lims.plot(x='time_nd',y='h1',linestyle='None',marker='s',color=colors[i],ax=ax,label=label,legend=False)
        label='{}_zenc'.format(L0_legend[case])
        df_lims.plot(x='time_nd',y='zenc',linestyle='-',color=colors[i],ax=ax,label=label,legend=False)
    ax.legend(numpoints=1,bbox_to_anchor=(1.,0.9),bbox_transform=fig.transFigure,prop={'size':11})
    ax.set_ylim([0,1800])
    ax.set_ylabel('height (m)')
    ax.set_title('h1 vs. zenc for all profiles')
fig.savefig('h1_time.png')

    
    
plt.show()

