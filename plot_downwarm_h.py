import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
import pandas as pd

from nchap_class import *
from nchap_class import For_Plots
from matplotlib import rcParams
rcParams.update({'font.size': 10})



"""
   for plotting the downward warm moving air at h
"""

def Main_Fun(rundate, gamma, flux_s):
     
     #output times
     dump_time_list, Times = Make_Timelists(1, 1800, 28800)
     Times = np.array(Times)  

     #class for pulling data files
     files = For_Plots(rundate)

     #Create lists of variable lists
     theta_file_list = [files.get_file(dump_time, "theta_bar") for dump_time in dump_time_list]

     press_file_list = [files.get_file(dump_time, "press") for dump_time in dump_time_list]
    
     flux_quads_file_list = [files.get_file(dump_time, "flux_quads_test") for dump_time in dump_time_list] #"flux_quads_theta1":[upwarm_bar, downwarm_bar, upcold_bar, downcold_bar, wvelthetapert_bar] 
     
     height_file = files.get_file("0000000600", "heights")

     #Get heights, scaling parameters
     AvProfVars = files.AvProfVars()
     rinovals = files.rinovals()
     gm_vars = files.gm_vars()
     #print AvProfVars.shape, rinovals.shape
     upwarm_h0 = []
     downwarm_h0=[]
     upcold_h0=[]
     downcold_h0=[]
     flux_h0=[]
     scaled_times=[]
     #loop over text files files
     for i in range(len(flux_quads_file_list)):
         #print i, theta_file_list[i]
         if rundate=="Nov302013": #different output frequency
             j = (i+1)*2 - 1
         else:
             j = (i+1)*3 - 1
         #print i, j
         height = np.genfromtxt(height_file)
         theta = np.genfromtxt(theta_file_list[i])
         #print theta.shape
         press = np.genfromtxt(press_file_list[i])
         rhow = nc.calc_rhow(press, height, theta[0])
         h = AvProfVars[j, 4]
         h0=AvProfVars[j, 4]
         h1 = AvProfVars[j, 2]
         scaled_time=np.divide(gm_vars[j, 3], gm_vars[j, 0])
         deltah = h1-h
         deltatheta = gamma*deltah
         thetastar = rinovals[j, 9]
         wstar = rinovals[j, 2]
         h0_lev = np.where(height == h0)
         flux_quads = np.genfromtxt(flux_quads_file_list[i])
         flux_s1 = 1.0*flux_s/(rhow[0]*1004)
         #flux_quads: 
         upwarm = flux_quads[h0_lev, 0][0][0]
         downwarm = flux_quads[h0_lev, 1][0][0]
         upcold = flux_quads[h0_lev, 2][0][0]
         downcold = flux_quads[h0_lev, 3][0][0]
         flux = flux_quads[h0_lev, 4][0][0]
         #print flux_s1, flux_s
         upwarm_h0.append(1.0*upwarm/(flux_s1))# // (thetastar) / /(gamma*deltah)/(0.2*thetastar)
         downwarm_h0.append(1.0*downwarm/(flux_s1))
         upcold_h0.append(1.0*upcold/(flux_s1))
         downcold_h0.append(1.0*downcold/(flux_s1))
         flux_h0.append(1.0*flux/(flux_s1))
         scaled_times.append(scaled_time)
     upwarm_h0 = np.array(upwarm_h0)
     upcold_h0 = np.array(upcold_h0)
     downcold_h0 = np.array(downcold_h0)
     downwarm_h0 = np.array(downwarm_h0)
     flux_h0 = np.array(flux_h0)

     if rundate=="Jan152014_1":
         scaled_times, flux_h0, upwarm_h0, downwarm_h0, upcold_h0, downcold_h0 = scaled_times[0:11], flux_h0[0:11], upwarm_h0[0:11], downwarm_h0[0:11], upcold_h0[0:11], downcold_h0[0:11] 
     
     return upwarm_h0, downwarm_h0, upcold_h0, downcold_h0, scaled_times, flux_h0    

run_list = [["Dec142013", .01, 100, '100/10', 'kv'], ["Nov302013", .005, 100, '100/5','ko'], ["Dec202013", .005, 60,'60/5','yo'], ["Dec252013", .0025, 60,'60/2.5','y*'], ["Jan152014_1", .005, 150, '150/5','ro'], ["Mar12014", .01, 60,'60/10','yv'], ["Mar52014", .01, 150,'150/10','rv']]

Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)
#Ax3.set_xlabel(r"$Time", fontsize=30)
Ax3.tick_params(axis="both", labelsize=20)
#Ax3.set_ylabel(r"$\frac{\overline{\theta^{\prime +}}_{h}}{\gamma ( h_{1}-h)} \ (where \ w^{\prime}<0)$", fontsize=30)
Ax3.set_ylabel(r"$\overline{w^{\prime}\theta^{\prime}}_{z_{f}}/\overline{w^{\prime}\theta^{\prime}}_{s}$", fontsize=30)
Ax3.set_xlabel(r"$z_{enc}/L_{0}$", fontsize=30)
#Ax3.set_ylabel(r"$\frac{\overline{w^{\prime-}_{h}}(where \ \theta^{\prime}>0) }{w^{*}}$ ", fontsize=30)

Ax3.text(5, .9, "(a)", fontsize=30)

#Ax3.set_ylabel( r"$\overline{w^{\prime -}}_{h}} \ (where \ \theta^{\prime}>0)", fontsize=30)

#Ax3.set_ylabel(r"$\overline{w^{\prime-}\theta^{\prime+}}_{h}$ (ms$^{-1}$K)", fontsize=30)

#Ax3.set_ylim(-.14, 0)
Ax3.set_ylim(-1.1, 1.1)
#Ax3.set_ylim(-.5, 0)
#Ax3.set_xlim(2, 8.2)

for run in run_list:
    #print run[0]
    upwarm_h0, downwarm_h0, upcold_h0, downcold_h0, scaled_times, flux_h0 = Main_Fun(run[0], run[1], run[2])
    downwarm_h0[0:1]=np.nan
    downcold_h0[0:1]=np.nan
    if run[0]=="Dec202013":
         downcold_h0[3:5]=np.nan
         downcold_h0[14]=np.nan
    Ax3.plot(scaled_times, 1.0*upwarm_h0, 'r+') #run[4], markersize=10
    Ax3.plot(scaled_times, 1.0*downwarm_h0, 'r--')# run[4]
    Ax3.plot(scaled_times, 1.0*upcold_h0, 'b+' )#run[4]
    Ax3.plot(scaled_times, 1.0*downcold_h0, 'b--')#run[4]
    Ax3.plot(scaled_times, 1.0*flux_h0, 'k-')#run[4]
#Ax3.legend(bbox_to_anchor=(1.49, 1.03), prop={'size':20}, numpoints = 1)
#box = Ax3.get_position()
#Ax3.set_position([box.x0, box.y0, box.width*1.33, box.height])
plt.tight_layout()
plt.show()























   
    
