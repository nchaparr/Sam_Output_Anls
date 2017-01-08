import h5py
from enum import IntEnum
from collections import OrderedDict as od
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import json

plt.style.use('ggplot')

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
# get all 6 heights and average theta from file produced
# by Flux_quads.py -- load thiese into a dictionary called var_dict
#
columns=list(Heights.__members__.keys())[:6]
var_dict = defaultdict(lambda: defaultdict(od))
h5_profiles='profiles_flux_quads_260.h5'
with h5py.File(h5_profiles,'r') as infile:
    case_list=list(infile.keys())
    for case in case_list:
        scales=infile[case]['scales'][...]
        scale_columns=json.loads(infile[case]['scales'].attrs['scale_columns'])
        avg_theta = infile[case]['ens_av_thetas'][...]
        avg_theta=np.mean(avg_theta,axis=(1,2))
        var_dict[case]['avg_theta']=avg_theta
        print('shape: ',avg_theta.shape)
        height = infile[case]['height'][...]
        hvals = infile[case]['hvals'][...]
        print('case: ',case)
        print('time index: ',infile[case]['time_list'].attrs['time_index'])
        time_index = int(infile[case]['time_list'].attrs['time_index'])
        var_dict[case]['wstar'] = scales[time_index,scale_columns['wstar']]
        var_dict[case]['thetastar']=scales[time_index,scale_columns['thetastar']]
        for col in columns:
            var_dict[case]['height_dict'][col]=hvals[time_index,int(Heights[col])]
        var_dict[case]['height']=height
        index_dict=od()
        for key in columns:
            height_value=var_dict[case]['height_dict'][key]
            index_dict[key]=np.searchsorted(var_dict[case]['height'],height_value,
                                            side='left')
        var_dict[case]['height_index']=index_dict

print('through with ',h5_profiles)        

quad_file='full_out_260.h5'
all_vars = ['wptp','wptn','wntn','wntp','wpert','thetapert']
quadrants = all_vars[:4]
#
# open the quadrant fluxes file produced by full_fluxes.py
# find the quadrant fluxes at each level  
#
with h5py.File(quad_file,'r') as flux_in:
    case_list=list(flux_in.keys())[:2]
    for case in case_list:
        print('adding heights for: ',case)
        top_lev=var_dict[case]['height_index']['zg1']
        bot_lev = var_dict[case]['height_index']['zg0']
        levs= list(range(bot_lev-1,top_lev+1))
        var_dict[case]['plot_heights']=var_dict[case]['height'][levs]
        var_dict[case]['zlevs']=levs
        for var in all_vars:
            for zlev in levs:
                szlev=str(zlev)
                var_dict[case][zlev][var]=flux_in[case][szlev][var][...]
                
print('through with ',quad_file)

                
def calc_stats(var):
    rms = float(np.sqrt(np.mean(var**2.)))
    std = float(np.std(var))
    count = len(var)
    return [rms,std,count]


op_dict=dict(wptp=(np.greater,np.greater),
             wptn=(np.greater,np.less),
             wntn=(np.less,np.less),
             wntp=(np.less,np.greater))

#
# for each level, recalculate the fluxes just to be safe using
# the logical comparisons from op_dict  -- we could also just use
# the values in var_dict
#
fieldnames=['Wrms','Wstd','Wcount','Trms','Tstd','Tcount']
record_list=[]
df_dict={}
for case in case_list:
    for lev in var_dict[case]['zlevs']:
        record=dict(case=case,lev=lev)
        wpert=var_dict[case][lev]['wpert']
        thetapert=var_dict[case][lev]['thetapert']
        record['tot_flux']=float(np.mean(wpert*thetapert))
        record['height'] = var_dict[case]['height'][lev]
        record['avg_theta']=var_dict[case]['avg_theta'][lev]
        for quad in quadrants: 
            op_w,op_t = op_dict[quad]
            hit=np.logical_and(op_w(wpert,0),op_t(thetapert,0))
            stat_vec=calc_stats(wpert[hit])
            Tout=calc_stats(thetapert[hit])
            stat_vec.extend(Tout)
            flux=float(np.mean(wpert[hit]*thetapert[hit]))
            group='{}_flux'.format(quad)
            record[group]=flux
            quad_names=['{}_{}'.format(quad,var) for var in fieldnames]
            stat_dict=dict(zip(quad_names,stat_vec))
            for field,value in stat_dict.items():
                record[field]=value
        record_list.append(record)
    df=pd.DataFrame(record_list)
    df_dict[case]=df
        
plt.close('all')

def plot2(ax1,df,the_key,height_dict,case_label='case'):
    """
     plot statistics from dataframe df, column the_key
     agains total flux
    """
    ax2 = ax1.twiny()
    ax1.plot('tot_flux','height','b-',data=df)
    ax2.plot(the_key,'height','r-',data=df,label=case_label)
    for t1 in ax1.get_xticklabels():
        t1.set_color('b')
    for t2 in ax2.get_xticklabels():
        t2.set_color('r')
    for key,hlev in height_dict.items():
        print('writing text: ',key,hlev)
        ax1.text(0.,hlev,key)
    ax1.set_xlabel('net flux (K m/s)')
    ax1.xaxis.label.set_color('b')
    ax2.xaxis.label.set_color('r')
    return ax1,ax2

quadrants=['wptp','wptn','wntp','wntn']


plot_dict=defaultdict(lambda: defaultdict(od))

plot_dict['W']['stats'] = ['Wrms','Wstd','Wcount']
plot_dict['T']['stats'] = ['Trms','Tstd','Tcount']
plot_dict['W']['labels'] = ['Wrms (K)','Wstd (K)','Wcount']
plot_dict['T']['labels'] = ['thetarms (K)','thestd (K)','thetacount']
plot_dict['T']['fig_title'] = 'theta'
plot_dict['W']['fig_title'] = 'wvel'

plt.close('all')
case_list=case_list[:3]
for quadrant in quadrants[:1]:
    print('working on ',case_list)
    for var in ['T','W']:
        stats=plot_dict[var]['stats']
        ax_labels=plot_dict[var]['labels']
        fig,axes=plt.subplots(2,2,figsize=(14,14))
        print('plotting: ',var,quadrant)
        flat_axes=axes.ravel()
        #
        # create variable name like wptp_Trms
        #
        plot_vars=['{}_{}'.format(quadrant,stat) for stat in stats]
        print('top of page: ',plot_vars,ax_labels)
        for count,plot_vals in enumerate(zip(plot_vars,ax_labels)):
            print('axes: ',plot_vals)
            plot_var,ax_label = plot_vals
            print('plotting variable: ',plot_var,ax_label)
            the_ax=flat_axes[count]
            casename=case_list[0]
            height_dict=var_dict[casename]['height_dict']
            print('first plot for ',casename)
            df_0=df_dict[casename]
            ax1,ax2 = plot2(the_ax,df_0,plot_var,height_dict,case_label=casename)
            # if len(case_list) > 1:
            #     fign,axn=plt.subplots(1,1)
            #     for the_case in case_list[1:]:
            #         df=df_dict[the_case]
            #         print('plotting line for ',the_case)
            #         #axn.plot(plot_var,'height',data=df,label=the_case)
            ax2.set_xlabel(ax_label)
            ax2.xaxis.label.set_color('r')
        casename=case_list[0]
        height_dict=var_dict[casename]['height_dict']
        the_ax=flat_axes[3]
        df_0=df_dict[casename]
        ax1,ax2 = plot2(the_ax,df_0,'avg_theta',height_dict)
        ax2.set_xlabel('avg_theta (K)')
        ax2.xaxis.label.set_color('r')
        fig_title='{} {} statistics'.format(quadrant,plot_dict[var]['fig_title'])
        fig_file='{}_fourplot_{}.png'.format(quadrant,plot_dict[var]['fig_title'])
        fig.suptitle(fig_title,size=20)
        fig.subplots_adjust(hspace=0.5)
        fig.savefig(fig_file)
        
plt.show()
