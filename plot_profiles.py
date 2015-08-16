import pandas as pd
import numpy as np
from collections import OrderedDict as od

def find_index(times,target):
    times=np.asarray(times)
    delta=np.diff(times)[0]
    index=np.where(np.abs(times - target) < delta/2.)[0][0]
    return index


h5file_new='good.h5'
#
# get the root series and dataframes
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
# find index for each case where non-dimensional time= 350 and 500
#
target_times=[350,500]
for nd_time in target_times:
    for case in case_list:
        Ntimes=df_dict[case,'AvProfLims']['Ntimes']
        times=df_dict[case,'AvProfLims']['times']
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
        real_time=times_dict[the_time,case]
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

flux_fig,flux_ax=plt.subplots(1,1)
grad_fig,grad_ax=plt.subplots(1,1)
for the_time in target_times:
    for case in compare_cases:
        real_time=times_dict[the_time,case]
        x=df_dict[case,'scaled_flux'][real_time]
        y=df_dict[case,'scaled_height'][real_time]
        label='{}_{}'.format(L0_legend[case],the_time)
        flux_ax.plot(x,y,label=label)
        x=df_dict[case,'scaled_dtheta'][real_time]
        grad_ax.plot(x,y,label=label)
flux_ax.set_title('flux profile at Ntime {} and {}'.format(*target_times))
grad_ax.set_title('gradient profile Ntime {} and {}'.format(*target_times))
[ax.set_ylim(0.65,0.85) for ax in [flux_ax,grad_ax]]
grad_ax.set_ylim(0.65,0.85)
flux_ax.set_xlim([-0.1,0.3])
grad_ax.set_xlim([-0.1,0.3])
flux_ax.legend(loc='best')
grad_ax.legend(loc='best')
filename='flux_comb.png'.format(the_time)
flux_fig.savefig(filename)
filename='grad_comb.png'.format(the_time)
grad_fig.savefig(filename)


n_colors=6
colors=plt.cm.rainbow(np.linspace(0,1,n_colors))
flux_fig,flux_ax=plt.subplots(1,1)
grad_fig,grad_ax=plt.subplots(1,1)
case_count=0
for the_time in target_times:
    for case in compare_cases:
        real_time=times_dict[the_time,case]
        print('real time',real_time)
        x=df_dict[case,'scaled_flux'][real_time]
        y=df_dict[case,'scaled_height'][real_time]
        label='{}_{}'.format(L0_legend[case],the_time)
        flux_ax.plot(x,y,label=label,color=colors[case_count])
        x=df_dict[case,'scaled_dtheta'][real_time]
        Ntimes=df_dict[case,'AvProfLims']['Ntimes']
        grad_ax.plot(x,y,label=label,color=colors[case_count])
        index=find_index(Ntimes,the_time)
        h0_height=df_dict[case,'AvProfLims']['h0'][index]
        zf0_height=df_dict[case,'AvProfLims']['zf0'][index]
        h_height=df_dict[case,'AvProfLims']['h'][index]
        grad_ax.axhline(h0_height/h_height,color=colors[case_count])
        flux_ax.axhline(zf0_height/h_height,color=colors[case_count])
        case_count+=1
flux_ax.set_title('zoomed flux profile at Ntime {} and {}'.format(*target_times))
grad_ax.set_title('zoomed gradient profile Ntime {} and {}'.format(*target_times))
[(ax.set_ylim(0.65,0.85),ax.set_ylabel('height/h')) for ax in [flux_ax,grad_ax]]
flux_ax.set_xlim([-0.1,0.3])
grad_ax.set_xlim([-0.1,0.3])
flux_ax.set_xlabel('flux/sfc_flux')
grad_ax.set_xlabel('gradient/gamma')
grad_ax.grid(True,which='both')
flux_ax.grid(True,which='both')
flux_ax.legend(loc='best')
grad_ax.legend(loc='best')
filename='flux_comb_zoom.png'.format(the_time)
flux_fig.savefig(filename)
filename='grad_comb_zoom.png'.format(the_time)
grad_fig.savefig(filename)


h0_fig,h0_ax=plt.subplots(1,1)
grad_count=0
flux_count=3
for case in compare_cases:
    h0_height=np.asarray(df_dict[case,'AvProfLims']['h0'])
    zf0_height=np.asarray(df_dict[case,'AvProfLims']['zf0'])
    h_height=np.asarray(df_dict[case,'AvProfLims']['h'])
    Ntimes=np.asarray(df_dict[case,'AvProfLims']['Ntimes'])
    label='{}_h0'.format(L0_legend[case])
    h0_ax.plot(Ntimes,h0_height/h_height,color=colors[grad_count],label=label)
    grad_count+=1
    label='{}_zf0'.format(L0_legend[case])
    h0_ax.plot(Ntimes,zf0_height/h_height,color=colors[flux_count],label=label)
    flux_count+=1
h0_ax.legend(loc='best')
h0_ax.set_xlim(300,600)
h0_ax.set_ylim(0.65,0.86)
h0_ax.grid(True)
h0_ax.set_xlabel('Ntime')
h0_ax.set_ylabel('h0/h and zf0/h')
h0_ax.set_title('scaled h0 and zf0 vs time for L0=18, 23, 29')
h0_fig.savefig('scaled_h0_zf0.png')

wrms_fig,wrms_ax=plt.subplots(1,1)
case_count=0
for the_time in target_times:
    for case in compare_cases:
        real_time=times_dict[the_time,case]
        print('real time',real_time)
        x=df_dict[case,'scaled_wrms'][real_time]
        y=df_dict[case,'scaled_height'][real_time]
        label='{}_{}'.format(L0_legend[case],the_time)
        wrms_ax.plot(x,y,label=label,color=colors[case_count])
        case_count+=1
wrms_ax.set_title('wrms/wstar profle Ntime {} and {}'.format(*target_times))
[(ax.set_ylim(0.65,0.85),ax.set_ylabel('height/h')) for ax in [wrms_ax]]
wrms_ax.set_xlim([0.3,0.6])
wrms_ax.set_xlabel('wrms/wstar')
wrms_ax.grid(True,which='both')
wrms_ax.legend(loc='best')
filename='wrms_comb.png'.format(the_time)
wrms_fig.savefig(filename)


plt.show()



