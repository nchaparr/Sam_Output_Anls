#ML_Heights
from __future__ import division
from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib import cm
from matplotlib import ticker
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sys
import blmfigs
from blmfigs import nchap_fun as nc
from blmfigs.Make_Timelist import Make_Timelists
import warnings
warnings.simplefilter('ignore', np.RankWarning)

"""
    Pulls Mixed Layer height determined by fastfit.pyx and does a hitogram
              
"""

datadir="/tera/phil/nchaparr/python/Plotting"
outputdir='/tera/phil/paper_figs/output'
invdir="data/invrinos"

dump_time_list, time_hrs = Make_Timelists(1, 3600, 28800)

date = "Mar52014"

rinovals = np.genfromtxt('{}/{}/{}'.format(datadir,date,invdir))

ml_height_hist_vars=[]

heightbins="data/heights0000028800"
height = np.genfromtxt('{}/{}/{}'.format(datadir,date,heightbins))
for i in range(len(dump_time_list)):
     dump_time=dump_time_list[i]
     if i == 4:
          ML_Heights=[]
          case_range = 10
          for case in range(case_range):
               mixed_layer_file='{}/{}/data/mixed_layer_height_{:d}_{}'.format(
                    datadir,date,(case+1),dump_time)
               print('reading {}'.format(mixed_layer_file))
               data=np.genfromtxt(mixed_layer_file)
               ML_Heights.append(data)

ML_Heights=np.array(ML_Heights)
[zvals, yvals, xvals] = ML_Heights.shape
ML_Heights = np.reshape(ML_Heights, (zvals*yvals, xvals))
height_bin_vols = nc.Bin_Peaks(ML_Heights, height)

theFig = plt.figure()
theFig.clf()
theAx = theFig.add_subplot(111)
#theAx.set_title('Histogram of local Mixed Layer Heights from 1 Case at 5 hrs')
theAx.set_xlabel('z (m)')
theAx.set_ylabel('count')
theAx.bar(height, height_bin_vols)     
theAx.set_xlim(0, 1500) #TODO:need to test axis limits first
theAx.set_ylim(0, 60000)
plt.show()
