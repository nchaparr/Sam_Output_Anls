import h5py
from enum import IntEnum
from collections import OrderedDict as od
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import pandas as pd
import json


# nd_time=23
# nd_time=25
# nd_time=18
# nd_time=30
matplotlib.use('Agg')
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

def calc_stats(var):
    rms = float(np.sqrt(np.mean(var**2.)))
    std = float(np.std(var))
    count = len(var)
    return [rms,std,count]


op_dict=dict(wptp=(np.greater,np.greater),
             wptn=(np.greater,np.less),
             wntn=(np.less,np.less),
             wntp=(np.less,np.greater))

import argparse
linebreaks=argparse.RawTextHelpFormatter
descrip=__doc__.lstrip()
parser = argparse.ArgumentParser(formatter_class=linebreaks,description=descrip)
parser.add_argument("--nd_time",type=str, help="zenc/l0")
args=parser.parse_args()
nd_time = args.nd_time
print('nd_time is: ',args.nd_time)

#
# get all 6 heights and average theta from file produced
# by Flux_quads.py -- load thiese into a dictionary called var_dict
#
columns=list(Heights.__members__.keys())[:6]
var_dict = defaultdict(lambda: defaultdict(od))
h5_profiles='profiles_flux_quads_{}.h5'.format(args.nd_time)

with h5py.File(h5_profiles,'r') as infile:
    case_list=list(infile.keys())[0:2]
    #case_list=case_list[::-1]
    for case in case_list:
        case_tups=enumerate(case_list)
        print('processing   {}'.format(case))
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

quad_file='full_out_{}.h5'.format(nd_time)
#quad_file='full_out_260.h5'
all_vars = ['wptp','wptn','wntn','wntp','wpert','thetapert']
quadrants = all_vars[:4]
#
# open the quadrant fluxes file produced by full_fluxes.py
# find the quadrant fluxes apt each level  
#
#print('above quad read: keys are: ',var_dict[case].keys())
with h5py.File(quad_file,'r') as flux_in:
    for case in case_list:
        case_dict=json.loads(flux_in.attrs['case_dict'])
        #case_list=list(flux_in.keys())
        print('adding heights for: ',case)
        top_lev=var_dict[case]['height_index']['zg1']
        bot_lev = var_dict[case]['height_index']['zg0'] - 15
        if bot_lev < 1: bot_lev=1
        levs= list(range(bot_lev-1,top_lev+1))
        print('retrieving these levs: ',levs)
        var_dict[case]['plot_heights']=var_dict[case]['height'][levs]
        var_dict[case]['zlevs']=levs
        for zlev in levs:
            var_dict[case][zlev]=defaultdict(lambda: defaultdict(od))
            for var in all_vars:
                szlev=str(zlev)
                var_dict[case][zlev]['quadfile'][var]=flux_in[case][szlev][var][...]
#print('through all reads, keys are: ',var_dict[case].keys())
            
# #
# # var_dict[case][zlev]['quadfile']['wptn']  is dimensional wptn flux
# #
# print('through with ',quad_file)

                

#
# for each level, recalculate the fluxes just to be safe using
# the logical comparisons from op_dict  -- we could also just use
# the values in var_dict
#
fieldnames=['Wrms','Wstd','Wcount','Trms','Tstd','Tcount']
record_list=[]
for case in case_list:
    try:
        zenc_l0=case_dict[case]['zenc_l0']
    except KeyError:
        print('exception at, zenc_l0 key missing for ',case)
        zenc_l0=np.sqrt(2*case_dict[case]['nd_time'])
        case_dict[case]['zenc_l0']=zenc_l0
    case_dict[case]['zenc']=zenc_l0*case_dict[case]['L0']
    var_dict[case]['height_dict']['zenc']=case_dict[case]['zenc']
    wstar = var_dict[case]['wstar']
    thetastar=var_dict[case]['thetastar']
    zg=var_dict[case]['height_dict']['zg']
    zenc = case_dict[case]['zenc']
    for lev in var_dict[case]['zlevs']:
        record=dict(case=case,lev=lev)
        wpert=var_dict[case][lev]['quadfile']['wpert']
        thetapert=var_dict[case][lev]['quadfile']['thetapert']
        record['tot_flux']=float(np.mean(wpert*thetapert))
        record['height'] = var_dict[case]['height'][lev]
        record['zg_nd_height']=var_dict[case]['height'][lev]/zg
        record['zenc_nd_height']=var_dict[case]['height'][lev]/zenc
        record['avg_theta']=var_dict[case]['avg_theta'][lev]
        for quad in quadrants: 
            op_w,op_t = op_dict[quad]
            hit=np.logical_and(op_w(wpert,0),op_t(thetapert,0))
            stat_vec=calc_stats(wpert[hit])
            Tout=calc_stats(thetapert[hit])
            var_dict[case][lev][quad]['wpert']=wpert[hit]
            var_dict[case][lev][quad]['thetapert']=thetapert[hit]
            stat_vec.extend(Tout)
            flux=float(np.mean(wpert[hit]*thetapert[hit]))
            group='{}_flux'.format(quad)
            record[group]=flux
            quad_names=['{}_{}'.format(quad,var) for var in fieldnames]
            stat_dict=dict(zip(quad_names,stat_vec))
            for field,value in stat_dict.items():
                record[field]=value
        record_list.append(record)
df_all=pd.DataFrame(record_list)

def get_perts(df_all,var_dict=None,bound=None,quad=None,case=None):
    df_0=df_all.loc[df_all['case'] == case]
    print('working with {} {}'.format(case,len(df_0)))
    target=var_dict[case]['height_dict'][bound]/var_dict[case]['height_dict']['zenc']
    zlevs=df_0['zenc_nd_height'].values
    index=int(np.searchsorted(zlevs,target))
    print('found: {} {} index: {} {}'.format(case,bound,index,len(df_0['lev'])))
    height_index=df_0['lev'].values[index]
    wpert=var_dict[case][height_index][quad]['wpert']
    thetapert=var_dict[case][height_index][quad]['thetapert']
    return wpert,thetapert

plt.close('all')
for case in case_list:
    fig,axes=plt.subplots(2,2,figsize=(13,13))
    axes = axes.ravel()
    quad='wptp'
    count=0
    for bound in ['zf0','zg0']:
        wpert,thetapert = get_perts(df_all,bound=bound,quad=quad,case=case,var_dict=var_dict)
        axes[count].hist(wpert,bins=30)
        axes[count].set(xlabel='wpert (m/s)',title='wpert {}'.format(bound))
        axes[count+1].hist(thetapert,bins=30)
        axes[count+1].set(xlabel='thetapert (K)',title='thetapert {}'.format(bound))
        count+=2
    
    title='{}  {}  time: {}'.format(quad,case,nd_time)
    fig.suptitle(title,fontsize=14)
    filename='hists_{}_{}.png'.format(case,nd_time)
    fig.savefig(filename)
    

