import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc

from nchap_class import *
from nchap_class import For_Plots
from matplotlib import rcParams
rcParams.update({'font.size': 10})



"""
   for plotting the downward warm moving air at h
"""

def Main_Fun(rundate, gamma, flux_s, the_label, the_legend):
     
     #output times
     dump_time_list, Times = Make_Timelists(1, 1800, 28800)
     Times = np.array(Times)  

     #class for pulling data files
     files = For_Plots(rundate)

     #Create lists of variable lists
     theta_file_list = [files.get_file(dump_time, "theta_bar") for dump_time in dump_time_list]

     press_file_list = [files.get_file(dump_time, "press") for dump_time in dump_time_list]
    
     flux_quads_file_list = [files.get_file(dump_time, "flux_quads_theta") for dump_time in dump_time_list]
     
     height_file = files.get_file("0000000600", "heights")

     AvProfVars = files.AvProfVars()
     rinovals = files.rinovals()
     print AvProfVars.shape, rinovals.shape
     downwarm_h = []
     #loop over text files files
     for i in range(len(flux_quads_file_list)):
         #print i, theta_file_list[i]
         if rundate=="Nov302013":
             j = (i+1)*2 - 1
         else:
             j = (i+1)*3 - 1
         print i, j
         height = np.genfromtxt(height_file)
         theta = np.genfromtxt(theta_file_list[i])
         #print theta.shape
         press = np.genfromtxt(press_file_list[i])
         rhow = nc.calc_rhow(press, height, theta[0])
         h = AvProfVars[j, 1]
         h0=AvProfVars[j, 0]
         h1 = AvProfVars[j, 2]
         deltah = h1-h
         deltatheta = gamma*deltah
         thetastar = rinovals[j, 9]
         wstar = rinovals[j, 2]
         h_lev = np.where(height == h)
         flux_quads = np.genfromtxt(flux_quads_file_list[i])
         flux_s1 = 1.0*flux_s/(rhow[0]*1004)         
         downwarm = flux_quads[h_lev, 1][0][0]
         print flux_s1, flux_s
         downwarm_h.append(1.0*downwarm/thetastar)

     downwarm_h = np.array(downwarm_h)    
     if rundate=="Jan152014_1":
         Times, downwarm_h = Times[0:11], downwarm_h[0:11] 
     Ax3.plot(Times, 1.0*downwarm_h, the_legend, label = the_label, markersize=12)

run_list = [["Dec142013", .01, 100, '100/10', 'kv'], ["Nov302013", .005, 100, '100/5','ko'], ["Dec202013", .005, 60,'60/5','yo'], ["Dec252013", .0025, 60,'60/2.5','y*'], ["Jan152014_1", .005, 150, '150/5','ro'], ["Mar12014", .01, 60,'60/10','yv'], ["Mar52014", .01, 150,'150/10','rv']]

Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)
Ax3.set_xlabel(r"$Time \ (hrs)$", fontsize=30)
Ax3.tick_params(axis="both", labelsize=20)
#Ax3.set_ylabel(r"$\frac{\overline{\theta^{\prime +}}_{h}}{\gamma ( h_{1}-h)} \ (where \ w^{\prime}<0)$", fontsize=30)
#Ax3.set_ylabel(r"$\frac{ \overline{w^{\prime-}\theta^{\prime+}}_{h}}{\overline{w^{\prime}\theta^{\prime}}_{s}}$", fontsize=30)
Ax3.set_ylabel(r"$\frac{\overline{\theta^{\prime +}}_{h} (where \ w^{\prime}<0) }{\theta^{*}}$", fontsize=30)
#Ax3.set_ylabel(r"$\frac{\overline{w^{\prime-}}(where \ \theta^{\prime}>0) }{w^{*}}$ ", fontsize=30)

#Ax3.set_ylabel(r"$\overline{\theta^{\prime+}}(where \ w^{\prime}<0) \ K$ ", fontsize=30)

#Ax3.set_ylabel( r"$\overline{w^{\prime -}}_{h}} \ (where \ \theta^{\prime}>0)", fontsize=30)

#Ax3.set_ylabel(r"$\overline{w^{\prime-}\theta^{\prime+}}_{h} \ (ms^{-1}K)$", fontsize=30)

Ax3.set_ylim(0, 2.2)
Ax3.set_xlim(2, 8.2)
for run in run_list:
    print run[0]
    Main_Fun(run[0], run[1], run[2], run[3], run[4])

Ax3.legend(loc = 'upper left', prop={'size':20}, numpoints = 1)
plt.tight_layout()
plt.show()


    
    
