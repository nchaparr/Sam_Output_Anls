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
#from matplotlib.patches import Patch
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from nchap_class import *
from nchap_class import For_Plots
from Make_Timelist import *
import warnings
warnings.simplefilter('ignore', np.RankWarning)
#import pywt



"""
    Pulls Mixed Layer height determined by fastfit.pyx and does a hitogram
              
"""


dump_time_list, time_hrs = Make_Timelists(1, 3600, 28800)
date_list = ["Dec252013"] #,, "Nov302013", "Jan152014_1",  
width_list = [.005] #, 2.5, 5, 10, .005, .005   
alpha_list = [.25] #, .5, 1   
label_list=[r"$\overline{w^{\prime}\theta^{\prime}}_{s} = 60$ Wm$^{-2}$"] #, r"$\overline{w^{\prime}\theta^{\prime}}_{s} = 100$ Wm$^{-2}$", r"$\overline{w^{\prime}\theta^{\prime}}_{s} = 150$ Wm$^{-2}$" 
color_list = ['0'] #,'k', '.5',  
theFig = plt.figure()
theFig.clf()

theAx = theFig.add_subplot(111)
theAx.set_title(r"$\gamma = 2.5$ (Kkm$^{-1}$)", fontsize=30)
theAx.set_ylabel(r'$P( \frac{h^{l}_{0}}{h} )$', fontsize=30)
#theAx.set_ylabel(r'$Count$', fontsize=30)
theAx.set_xlabel(r'$\frac{h^{l}_{0}}{h}$', fontsize=32)
#theAx.set_xlabel(r'$h^{l}_{0}$ (m)', fontsize=32)
theAx.set_xticks([.6, .8, 1, 1.2])
theAx.set_yticks([0, 0.05, 0.1])
#theAx.set_xticks([500, 1000, 1500])
#theAx.set_yticks([0, 10000, 20000])
#theAx.set_xlim(400, 1500) #TODO:need to test axis limts first
#theAx.set_ylim(0, 25000)

theAx.set_xlim(0.6, 1.2) #TODO:need to test axis limts first
theAx.set_ylim(0, .1)

#theAx.text(1060, 10000, r"$\overline{w^{\prime}\theta^{\prime}}_{s} = 150 \ Wm^{-2}$", fontsize=20)
#theAx.text(900, 12000, r"$\overline{w^{\prime}\theta^{\prime}}_{s} = 100 \ Wm^{-2}$", fontsize=20)
#theAx.text(0.9, .04, r"$\overline{w^{\prime}\theta^{\prime}}_{s} = 60 \ Wm^{-2}$", fontsize=20)
theAx.tick_params(axis='both', which='major', labelsize=30)
#rinovals = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/invrinos")
for j in range(len(date_list)):
    date=date_list[j]
    Width=width_list[j]
    Color = color_list[j]
    Alpha = alpha_list[j]
    Label = label_list[j]
    ml_height_hist_vars=[]
    points = For_Plots(date)
    rinovals = points.rinovals()
    AvProfVars = points.AvProfVars()
         
    for i in range(len(dump_time_list)):
        dump_time=dump_time_list[i]
        if i == 4:
            ML_Heights=[]
            case_range = 10
            for case in range(case_range):
                height = np.genfromtxt("/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/heights0000028800")
                ML_Heights.append(np.genfromtxt("/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/mixed_layer_height_"+ str(case+1) + "_" + dump_time))

            ML_Heights=np.array(ML_Heights)
            [zvals, yvals, xvals] = ML_Heights.shape
            ML_Heights = np.reshape(ML_Heights, (zvals*yvals*xvals,))
            v_max, v_min, mean, var = np.amax(ML_Heights), np.amin(ML_Heights), np.mean(ML_Heights), np.var(ML_Heights)
          #print 'max min std', v_max, v_min, mean, var
            if date=="Nov302013":
                rinovals_index = (i+1)*4-1
            else:
                rinovals_index = (i+1)*6-1
          #print rinovals_index
          #ml_height_hist_vars.append([rinovals[rinovals_index,1], rinovals[rinovals_index,3], v_max, v_min, mean, var])
          #n, bins, patches = theAx.hist(tracer_peaks, bins=20)
            ML_Heights = np.reshape(ML_Heights, (zvals*yvals, xvals))
            height_bin_vols = nc.Bin_Peaks(ML_Heights, height)

          #set up plot
            h = AvProfVars[rinovals_index, 1]
          #theAx.set_title('Histogram of local Mixed Layer Heights from 1 Case at 5 hrs')
            theAx.bar(1.0*height/h, 1.0*height_bin_vols/(zvals*yvals*xvals), width = Width, alpha = Alpha, color=Color, label=Label)#h  /     
#theAx.legend(prop={'size':24.5}) #bbox_to_anchor=(1.45, 1),, borderaxespad=0.0
#box = theAx.get_position()
#theAx.set_position([box.x0, box.y0, box.width*1.2, box.height])
theFig.tight_layout()          
theFig.show()          
#theFig.savefig("/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/ML_Height_hist.png")

#np.savetxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/ml_height_hist_vars", np.array(ml_height_hist_vars), delimiter=' ')

    
    
