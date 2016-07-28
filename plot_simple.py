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
from nchap_fun import Get_CBLHeights
import numpy as np
import pdb
import pandas as pd
import errno, sys


if __name__ == "__main__":

    linebreaks=argparse.RawTextHelpFormatter
    descrip = __doc__.lstrip()
    parser = argparse.ArgumentParser(formatter_class=linebreaks,description=descrip)
    parser.add_argument('h5_file', type=str,help='h5 file to plot')
    args=parser.parse_args()


    with pd.HDFStore(args.h5_file,'r') as store:
        df_old = store.get('new/limits')
        df_new = store.get('old/limits')

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
