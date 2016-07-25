"""
 plot profiles form an h5 file produced by get_thetas.py

 example:

 python plot_simple.py thetaprofs_Dec142013.h5
"""

from matplotlib import pyplot as plt
import h5py
import json
from collections import OrderedDict
from matplotlib.colors import Normalize
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import argparse
from nchap_fun import Get_CBLHeights
import numpy as np
import pdb

linebreaks=argparse.RawTextHelpFormatter
descrip = __doc__.lstrip()
parser = argparse.ArgumentParser(formatter_class=linebreaks,description=descrip)
parser.add_argument('h5_file', type=str,help='h5 file to plot')
args=parser.parse_args()

with h5py.File(args.h5_file,'r') as f:
    case_dict=json.loads(f.attrs['case_dict'],object_pairs_hook=OrderedDict)
    thetas = f['thetas'][...]
    heat_flux = f['fluxes'][...]
    height = f['height'][...]
    press = f['press'][...]*0.01  #convert to pascals

rundate = case_dict['name']
gammas = case_dict['gammas']
flux_s = case_dict['fluxes']

if rundate == "Jan152014_1":
    top_index = np.where(abs(2000 - height) < 26.)[0][0] #may need to be higher (e.g. for 60/2.5)
else:
    top_index = np.where(abs(1700 - height) < 26.)[0][0] #may need to be higher (e.g. for 60/2.5)

    
the_times = case_dict['float_hours']

ntimes,nheights = thetas.shape

pdb.set_trace()
for timestep in range(ntimes):
    theta_prof = thetas[timestep,:]
    flux_prof = heat_flux[timestep,:]
    #[elbot_dthetadz, h, eltop_dthetadz, elbot_flux ,h_flux  ,eltop_flux, deltatheta, mltheta, z1_GM]= \
    ez_heights =  Get_CBLHeights(height, press, theta_prof, flux_prof, gammas, flux_s, top_index, 'new')


plt.close('all')
fig, ax = plt.subplots(1,1)

colors=sns.color_palette('coolwarm')
pal=LinearSegmentedColormap.from_list('test',colors)
pal.set_bad('0.75') #75% grey
pal.set_over('r')
pal.set_under('k')
vmin= 299.
vmax= 310.
the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
image = ax.pcolormesh(the_times,height,thetas.T,cmap=pal,norm=the_norm)
cax=fig.colorbar(image)
ax.set(ylim=(0,1000))

ntimes,nheights = thetas.shape


fig, ax = plt.subplots(1,1)
for row in range(ntimes):
    theta = thetas[row,:]
    ax.plot(theta,height)
ax.set(xlim=(300,308),ylim=(0,1000))


fig, ax = plt.subplots(1,1)
ntimes,nheights = heat_flux.shape
for row in range(ntimes):
    flux = heat_flux[row,:]
    ax.plot(flux,height)
ax.set(ylim=(0,1000))

plt.show()
