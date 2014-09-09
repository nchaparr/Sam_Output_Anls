import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from nchap_class import *
from matplotlib import rcParams
rcParams.update({'font.size': 10})



"""
   for plotting the downward warm moving air at h
"""

def Main_Fun(rundate, gamma, flux_s):
     
     #output times
     dump_time_list, Times = Make_Timelists(1, 600, 28800)
     Times = np.array(Times)  

     #class for pulling data files
     files = For_Plots(rundate)

     #Create lists of variable lists
     
     flux_quads_file_list = [files.get_file(dump_time, "flux_quads") for dump_time in dump_time_list]
     
     height_file = files.get_file("0000000600", "heights")

     AvProfVars = files.AvProfVars()
     downwarm_h = []
     #loop over text files files
     for i in range(len(flux_quad_file_list)):
         #print i, theta_file_list[i]
         j = (i+1)*3 - 1
         print i, j
         theta = np.genfromtxt(theta_file_list[i])
         #print theta.shape
         height = np.genfromtxt(height_file)
         h = AvProfVars[j, 1]
         h_lev = np.where(height == h)
         flux_quads = np.genfromtxt(flux_quad_file_list[i])
         downwarm = flux_quads[h_lev]

run_list = [["Nov302013", .005, 100], ["Dec142013", .01, 100], ["Dec202013", .005, 60], ["Dec252013", .0025, 60], ["Jan152014_1", .005, 150], ["Mar12014", .01, 60], ["Mar52014", .01, 150]]

for run in run_list:
    print run
    Main_Fun(run[0], run[1], run[2])


    
    
