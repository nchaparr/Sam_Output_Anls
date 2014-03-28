from __future__ import division
from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib import cm
from matplotlib import ticker
import matplotlib.pyplot as plt
#import site
#site.addsitedir('/tera/phil/nchaparr/SAM2/sam_main/python')
#from Percentiles import *
from matplotlib.patches import Patch
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from Make_Timelist import *
import warnings
warnings.simplefilter('ignore', np.RankWarning)
#import pywt



"""
    Pulls Mixed Layer height determined by fastfit.pyx and does a hitogram
              
"""


dump_time_list, time_hrs = Make_Timelists(1, 3600, 28800)

date = "Nov302013"

rinovals = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/invrinos")

ml_height_hist_vars=[]

for i in range(len(dump_time_list)):
     dump_time=dump_time_list[i]
     if i == 4:
          ML_Heights=[]
          case_range = 10
          for case in range(case_range):
               height = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/heights0000028800")
               ML_Heights.append(np.genfromtxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/mixed_layer_height_"+ str(case+1) + "_" + dump_time))

          ML_Heights=np.array(ML_Heights)
          [zvals, yvals, xvals] = ML_Heights.shape
          ML_Heights = np.reshape(ML_Heights, (zvals*yvals*xvals,))
          v_max, v_min, mean, var = np.amax(ML_Heights), np.amin(ML_Heights), np.mean(ML_Heights), np.var(ML_Heights)
          print 'max min std', v_max, v_min, mean, var
          rinovals_index = (i+1)*6-1
          ml_height_hist_vars.append([rinovals[rinovals_index,1], rinovals[rinovals_index,3], v_max, v_min, mean, var])
          #n, bins, patches = theAx.hist(tracer_peaks, bins=20)
          ML_Heights = np.reshape(ML_Heights, (zvals*yvals, xvals))
          height_bin_vols = nc.Bin_Peaks(ML_Heights, height)
          #set up plot

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
          theFig.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/ML_Height_hist"+str(time_hrs[i])+".png")

np.savetxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/ml_height_hist_vars", np.array(ml_height_hist_vars), delimiter=' ')

    
    
