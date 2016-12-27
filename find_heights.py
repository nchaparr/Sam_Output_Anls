import h5py
from enum import IntEnum
from collections import OrderedDict as od
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt
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
results_dict = defaultdict(lambda: defaultdict(od))
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
    rms = np.sqrt(np.mean(var**2.))
    std = np.std(var)
    count = len(var)
    return rms,std,count
    
    
fieldnames=['wrms','wstd','wcount','Trms','Tstd','Tcount']
for case in case_list:
    for lev in var_dict[case]['zlevs']:
        wpert=var_dict[case][lev]['wpert']
        thetapert=var_dict[case][lev]['thetapert']
        results_dict[case][lev]['flux'] = np.mean(wpert*thetapert)
        for quad in quadrants: 
            op_w,op_t = op_dict[key]
            hit=np.logical_and(op_w(wpert,0),op_t(thetapert,0))
            wout=calc_stats(wpert[hit])
            Tout=calc_stats(thetapert[hit])
            stat_vec=wout.extend(Tout)
            stat_dict=dict(zip(fieldnames,stat_vec))
            results_dict[case][lev][quad].update(stat_dict)

#for case in case_list:
#             for zlev in levs:
#                 print('zlev: ',zlev)
#                 wpert=flux_in[case][str(zlev)]['wpert'][...]
#                 thetapert=flux_in[case][str(zlev)]['thetapert'][...]
#                 uw_hit = ('wptp',np.logical_and(thetapert >0, wpert >0))
#                 uc_hit = ('wptn',np.logical_and(thetapert <0, wpert >0))
#                 dc_hit = ('dc',np.logical_and(thetapert <0, wpert <0))
#                 dw_hit = ('dw',np.logical_and(thetapert >0, wpert <0))
                
                
#                 theta_uw=thetapert[uw_hit]
#                 tot_flux=(wpert*thetapert).mean()
#                 for var in ['mean_flux','theta_rms','theta_std','theta_count']:
#                     theta_out = np.sqrt(np.mean(theta_uw**2.))
#                     theta_std = np.std(theta_uw)
#                     theta_count = len(theta_uw)
#                     var_dict['mean_flux'].append(flux)
#                     var_dict['theta_rms'].append(theta_out)
#                     var_dict['theta_std'].append(theta_std)
#                     var_dict['theta_count'].append(theta_count)
        
# plt.close()

# def plot2(ax1,var_dict,key1,key2,zlevs,height_dict):
#     prof1=var_dict[key1]
#     prof2=var_dict[key2]
#     ax2 = ax1.twiny()
#     ax1.plot(prof1,zlevs,'b-')
#     ax2.plot(prof2,zlevs,'r-')
#     for t1 in ax1.get_xticklabels():
#         t1.set_color('b')
#     for t2 in ax2.get_xticklabels():
#         t2.set_color('r')
#     for key,hlev in height_dict.items():
#         ax1.text(0.,hlev,key)
#     ax1.set_xlabel('net flux (K m/s)')
#     ax1.xaxis.label.set_color('b')
#     ax2.xaxis.label.set_color('r')
#     return ax1,ax2


# zlevs=height[levs]
# avg_theta_plot=avg_theta[levs]
# var_dict['avg_theta']=avg_theta
# fig,axes=plt.subplots(2,2,figsize=(14,14))
# ax1,ax2 = plot2(axes[0,0],var_dict,'mean_flux','theta_rms',zlevs,height_dict)
# ax2.set_xlabel('thetarms (K)')
# ax2.xaxis.label.set_color('r')

# ax1,ax2 = plot2(axes[0,1],var_dict,'mean_flux','theta_std',zlevs,height_dict)
# ax2.set_xlabel('thetstd (K)')

# ax1,ax2 = plot2(axes[1,0],var_dict,'mean_flux','theta_count',zlevs,height_dict)
# ax2.set_xlabel('thetacount')

# ax1,ax2 = plot2(axes[1,1],var_dict,'mean_flux','avg_theta',zlevs,height_dict)
# ax2.set_xlabel('avg_theta')


# fig.subplots_adjust(hspace=0.5)
# fig.savefig('fourplot.png')
        
# plt.show()
