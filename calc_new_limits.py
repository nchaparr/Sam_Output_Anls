"""
 plot profiles form an h5 file produced by get_thetas.py

 example:

 python plot_simple.py thetaprofs_Mar12014.h5

 or

 python plot_simple.py thetaprofs_Mar12014.h5 --new
"""

from matplotlib import pyplot as plt
import h5py, os
import json
from collections import OrderedDict
from matplotlib.colors import Normalize
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import argparse
from nchap_fun import Get_CBLHeights, gm_vars
import numpy as np
import pdb
import pandas as pd
import errno, sys

def find_heights(h5_file,old_new):
    if os.path.exists(args.h5_file):
        with h5py.File(args.h5_file,'r') as f:
            case_dict=json.loads(f.attrs['case_dict'],object_pairs_hook=OrderedDict)
            thetas = f['thetas'][...]
            heat_flux = f['mean_fluxes'][...]
            height = f['height'][...]
            press = f['press'][...]*0.01  #convert to pascals
    else:
        print("{} isn't one of {}".format(args.h5_file,case_list))
        sys.exit(1)

    rundate = case_dict['name']
    gammas = case_dict['gammas']*1.e-3  #K/m
    flux_s = case_dict['fluxes']  #W/m^2
    the_times = np.array(case_dict['float_hours'])
    L0, N, B0, zenc = gm_vars(the_times, flux_s,gammas)
    nd_time = the_times*3600.*N
    gm_time = zenc/L0
    t_columns = ['t_hours','ndtime','gmtime']
    df_times = pd.DataFrame.from_records(zip(the_times,nd_time,gm_time),columns = t_columns)
    
    if rundate == "Jan152014_1":
        top_index = np.where(abs(2000 - height) < 26.)[0][0] #may need to be higher (e.g. for 60/2.5)
    else:
        top_index = np.where(abs(1700 - height) < 26.)[0][0] #may need to be higher (e.g. for 60/2.5)



    ntimes,nheights = thetas.shape

    #pdb.set_trace()
    height_vec = []
    for timestep in range(ntimes):
        theta_prof = thetas[timestep,:]
        flux_prof = heat_flux[timestep,:]
        #[elbot_dthetadz, h, eltop_dthetadz, elbot_flux ,h_flux  ,eltop_flux, deltatheta, mltheta, z1_GM]= \
        ez_heights =  Get_CBLHeights(height, press, theta_prof, flux_prof, gammas, flux_s, top_index, old_new)
        height_vec.append(ez_heights)

    
    columns=['zg0','zg','zg1','zf0','zf','zf1','deltatheta','mltheta','z1_GM']        
    df=pd.DataFrame.from_records(height_vec,columns=columns)
    comb_df = pd.merge(df_times,df,left_index=True, right_index=True)
    pdb.set_trace()
    return comb_df

if __name__ == "__main__":

    linebreaks=argparse.RawTextHelpFormatter
    descrip = __doc__.lstrip()
    parser = argparse.ArgumentParser(formatter_class=linebreaks,description=descrip)
    parser.add_argument('h5_file', type=str,help='h5 file to plot')
    parser.add_argument('--new',dest='new',action='store_true',help='use new scaling')
    args=parser.parse_args()

    old_new = 'old'
    if args.new:
        old_new = 'new'

    case_list = ['Dec252013', 'Mar12014', 'Jan152014_1', 'Nov302013', 'Dec142013', 'Mar52014', 'Dec202013']

    df = find_heights(args.h5_file,old_new)

    if args.new:
        with pd.HDFStore(args.h5_file,'a') as store:
            store.put('new/limits',df)
    else:
        with pd.HDFStore(args.h5_file,'a') as store:
            store.put('old/limits',df)

    print(df)

# plt.close('all')
# fig, ax = plt.subplots(1,1)
#the_times = case_dict['float_hours']

# colors=sns.color_palette('coolwarm')
# pal=LinearSegmentedColormap.from_list('test',colors)
# pal.set_bad('0.75') #75% grey
# pal.set_over('r')
# pal.set_under('k')
# vmin= 299.
# vmax= 310.
# the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
# image = ax.pcolormesh(the_times,height,thetas.T,cmap=pal,norm=the_norm)
# cax=fig.colorbar(image)
# ax.set(ylim=(0,1000))

# ntimes,nheights = thetas.shape


# fig, ax = plt.subplots(1,1)
# for row in range(ntimes):
#     theta = thetas[row,:]
#     ax.plot(theta,height)
# ax.set(xlim=(300,308),ylim=(0,1000))


# fig, ax = plt.subplots(1,1)
# ntimes,nheights = heat_flux.shape
# for row in range(ntimes):
#     flux = heat_flux[row,:]
#     ax.plot(flux,height)
# ax.set(ylim=(0,1000))

# plt.show()
