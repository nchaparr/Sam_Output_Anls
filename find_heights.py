import h5py
from enum import IntEnum
from collections import OrderedDict as od
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
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
# get all 6 heigths
#
columns=list(Heights.__members__.keys())[:6]
var_dict = defaultdict(lambda: defaultdict(od))
h5_profiles='profiles_flux_quads_260.h5'
with h5py.File(h5_profiles,'r') as infile:
    case_list=list(infile.keys())[:1]
    for case in case_list:
        avg_theta = infile[case]['ens_av_thetas'][...]
        avg_theta=np.mean(avg_theta,axis=(1,2))
        var_dict[case]['avg_theta']=avg_theta
        print('shape: ',avg_theta.shape)
        height = infile[case]['height'][...]
        hvals = infile[case]['hvals'][...]
        print('case: ',case)
        print('time index: ',infile[case]['time_list'].attrs['time_index'])
        time_index = int(infile[case]['time_list'].attrs['time_index'])
        for col in columns:
            var_dict[case][col]=hvals[time_index,int(Heights[col])]
        var_dict[case]['height']=height
        index_dict=od()
        for key in columns:
            height_value=var_dict[case][key]
            index_dict[key]=np.searchsorted(var_dict[case]['height'],height_value,
                                            side='left')
        var_dict[case]['height_index']=index_dict

prof_file='full_out_260.h5'
all_vars = ['wptp','wptn','wntn','wntp','wpert','thetapert']
quadrants = all_vars[:4]
op_dict=dict(wptp=(np.greater,np.greater),
             wptn=(np.greater,np.less),
             wntn=(np.less,np.less),
             wntp=(np.less,np.greater))
              
with h5py.File(prof_file,'r') as flux_in:
    case_list=list(flux_in.keys())[:1]
    for case in case_list:
        top_lev=var_dict[case]['height_index']['zg1']
        bot_lev = var_dict[case]['height_index']['zg0']
        levs= list(range(bot_lev-1,top_lev+1))
        var_dict[case]['plot_heights']=var_dict[case]['height'][levs]
        var_dict[case]['zlevs']=levs
        for var in all_vars:
            for zlev in levs:
                szlev=str(zlev)
                var_dict[case][zlev][var]=flux_in[case][szlev][var][...]

def calc_stats(var):
    rms = float(np.sqrt(np.mean(var**2.)))
    std = float(np.std(var))
    count = len(var)
    return [rms,std,count]
    
    
fieldnames=['wrms','wstd','wcount','Trms','Tstd','Tcount']
record_list=[]
for case in case_list:
    for lev in var_dict[case]['zlevs']:
        record=dict(case=case,lev=lev)
        wpert=var_dict[case][lev]['wpert']
        thetapert=var_dict[case][lev]['thetapert']
        record['tot_flux']=float(np.mean(wpert*thetapert))
        record['height'] = var_dict[case]['height'][lev]
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
print(df.head())            
        
plt.close()

def plot2(ax1,df,the_key,height_dict):
    ax2 = ax1.twiny()
    ax1.plot('tot_flux','height','b-',data=df)
    ax2.plot(the_key,'height','r-')
    for t1 in ax1.get_xticklabels():
        t1.set_color('b')
    for t2 in ax2.get_xticklabels():
        t2.set_color('r')
    for key,hlev in height_dict.items():
        ax1.text(0.,hlev,key)
    ax1.set_xlabel('net flux (K m/s)')
    ax1.xaxis.label.set_color('b')
    ax2.xaxis.label.set_color('r')
    return ax1,ax2


# zlevs=height[levs]
# avg_theta_plot=avg_theta[levs]
# var_dict['avg_theta']=avg_theta
# fig,axes=plt.subplots(2,2,figsize=(14,14))
# ax1,ax2 = plot2(axes[0,0],var_dict,'mean_flux','theta_rms',zlevs,height_dict)
# ax2.set_xlabel('thetarms (K)')
# ax2.xaxis.label.set_color('r')

# # ax1,ax2 = plot2(axes[0,1],var_dict,'mean_flux','theta_std',zlevs,height_dict)
# # ax2.set_xlabel('thetstd (K)')

# # ax1,ax2 = plot2(axes[1,0],var_dict,'mean_flux','theta_count',zlevs,height_dict)
# # ax2.set_xlabel('thetacount')

# # ax1,ax2 = plot2(axes[1,1],var_dict,'mean_flux','avg_theta',zlevs,height_dict)
# # ax2.set_xlabel('avg_theta')


# # fig.subplots_adjust(hspace=0.5)
# # fig.savefig('fourplot.png')
        
# # plt.show()
