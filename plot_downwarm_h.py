import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
#import pandas as pd

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
    
     flux_quads_file_list = [files.get_file(dump_time, "upwarm_rtmnsq_thetas") for dump_time in dump_time_list] #"flux_quads_theta1":[upwarm_bar, downwarm_bar, upcold_bar, downcold_bar, wvelthetapert_bar]
     flux_quads_file_list1 = [files.get_file(dump_time, "upwarm_rtmnsq_wvel") for dump_time in dump_time_list] #"flux_quads_theta1":[upwarm_bar, downwarm_bar, upcold_bar, downcold_bar, wvelthetapert_bar]
     
     height_file = files.get_file("0000000600", "heights")

     #Get heights, scaling parameters
     AvProfVars = files.AvProfVars()
     rinovals = files.rinovals()
     gm_vars = files.gm_vars()
     #print AvProfVars.shape, rinovals.shape
     upwarm_temp_h0 = []
     upwarm_wvel_h0 = []
     scaled_time = []
     #downwarm_h0=[]
     #upcold_h0=[]
     #downcold_h0=[]
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
         L0=gm_vars[j,0]
         zenc=gm_vars[j,3]
         h = AvProfVars[j, 4]
         h0=AvProfVars[j, 1]
         h1 = AvProfVars[j, 2]
         deltah = h1-h
         deltatheta = gamma*deltah
         thetastar = rinovals[j, 9]
         wstar = rinovals[j, 2]
         h0_lev = np.where(height==h0)
         #h0_lev = np.where(np.abs(h/2-height)<25)[0]
         flux_quads = np.genfromtxt(flux_quads_file_list[i])
         flux_quads1=np.genfromtxt(flux_quads_file_list1[i])
         flux_s1 = 1.0*flux_s/(rhow[0]*1004)
         #flux_quads: 
         upwarm_temp = flux_quads[h0_lev][0]
         upwarm_vel = flux_quads1[h0_lev][0]
         #downwarm = flux_quads[h0_lev, 1][0][0]
         #upcold = flux_quads[h0_lev, 2][0][0]
         #downcold = flux_quads[h0_lev, 3][0][0]
         #print flux_s1, flux_s
         upwarm_temp_h0.append(1.0*upwarm_temp/thetastar)# // (thetastar) / /(gamma*deltah)/(0.2*thetastar)
         upwarm_wvel_h0.append(1.0*upwarm_vel/wstar)# // (thetastar) / /(gamma*deltah)/(0.2*thetastar)
         scaled_time.append(1.0*zenc/L0)
         #downwarm_h0.append(1.0*downwarm/(flux_s1))
         #upcold_h0.append(1.0*upcold/(flux_s1))
         #downcold_h0.append(1.0*upcold/(flux_s1))

     upwarm_temp_h0 = np.array(upwarm_temp_h0)
     upwarm_wvel_h0 = np.array(upwarm_wvel_h0)
     #upcold_h0 = np.array(upcold_h0)
     #downcold_h0 = np.array(downcold_h0)
     #downwarm_h0 = np.array(downwarm_h0)

     if rundate=="Jan152014_1":
         scaled_time, upwarm_temp_h0, upwarm_wvel_h0 = scaled_time[0:11], upwarm_temp_h0[0:11], upwarm_wvel_h0[0:11] 
     
     return upwarm_temp_h0, upwarm_wvel_h0, scaled_time    

run_list = [["Dec142013", .01, 100, '100/10', 'kv'], ["Nov302013", .005, 100, '100/5','ko'], ["Dec202013", .005, 60,'60/5','yo'], ["Dec252013", .0025, 60,'60/2.5','y*'], ["Jan152014_1", .005, 150, '150/5','ro'], ["Mar12014", .01, 60,'60/10','yv'], ["Mar52014", .01, 150,'150/10','rv']]

Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)
Ax3.set_xlabel(r"$z_{enc}/L_{0}$", fontsize=30)
#Ax3.set_ylabel(r"$rms / \theta^{\prime}, w^{\prime}$", fontsize=30)
Ax3.set_ylabel(r"$\sqrt{\overline{w^{\prime 2}}}/w_{*}$", fontsize=30)
#Ax3.set_ylabel(r"$\sqrt{\overline{\theta^{\prime}}}$", fontsize=30)
Ax3.tick_params(axis="both", labelsize=20)
#Ax3.set_ylabel(r"$\sqrt{\overline{(w^{\prime})^{2}}}/w_*$", fontsize=30)
#Ax3.set_ylabel(r"$\frac{ \overline{w^{\prime-}\theta^{\prime+}}_{h}}{\overline{w^{\prime}\theta^{\prime}}_{s}}$", fontsize=30)
#Ax3.set_ylabel(r"$\frac{\overline{\theta^{\prime +}}_{h} (where \ w^{\prime}<0) }{\theta^{*}}$", fontsize=30)
#Ax3.set_ylabel(r"$\frac{\overline{w^{\prime-}_{h}}(where \ \theta^{\prime}>0) }{w^{*}}$ ", fontsize=30)
#Ax3.text(5, .01, r'(b)', fontsize=30)
#Ax3.set_ylim(0,2)
for run in run_list: 
     if run[0]=="Mar12014":
          upwarm_temp_h0, upwarm_wvel_h0, Times = Main_Fun(run[0], run[1], run[2])
          #Ax3.plot(Times, 1.0*upwarm_temp_h0, run[4], markersize=10)
          Ax3.plot(Times, 1.0*upwarm_wvel_h0, run[4], markersize=8)
     else:
          upwarm_temp_h0, upwarm_wvel_h0, Times = Main_Fun(run[0], run[1], run[2])
          #Ax3.plot(Times, 1.0*upwarm_temp_h0, run[4], markersize=10)
          Ax3.plot(Times, 1.0*upwarm_wvel_h0, run[4], markersize=8)
    #
    #Ax3.plot(Times, 1.0*downwarm_h0, run[4])
    #Ax3.plot(Times, 1.0*upcold_h0, 'b-')
    #Ax3.plot(Times, 1.0*downcold_h0, 'b-')
Ax3.text(5, .01, r'(b)', fontsize=30)
#Ax3.legend(loc="upper right", prop={'size':20}, numpoints = 1)
#box = Ax3.get_position()
#Ax3.set_position([box.x0, box.y0, box.width*1.33, box.height])
plt.tight_layout()
plt.show()


    
    
