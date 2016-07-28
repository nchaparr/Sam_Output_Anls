import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import Make_Timelists
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from nchap_class import For_Plots
from matplotlib import rcParams
import pdb
rcParams.update({'font.size': 10})

"""calculates temperature gradients (discrete) from txt files inturn from ensemble run 3D files
   gets levels where gradient exceeds zero, and where it resumes gamma, and zero crossings for fluxes
   and the maximal points withing the entrainment region.
   Gets Delta Theta and the Mixed Layer average Theta.
   Dumps them in a text file.
   Calcs and dumps rino, invrino, wstar
"""

def Main_Fun(rundate, gamma, flux_s):
     
     #output times
       
     if rundate=="Nov302013":
          dump_time_list, Times = Make_Timelists(1, 900, 28800)
          Times = np.array(Times)  
     else:
          dump_time_list, Times = Make_Timelists(1, 600, 28800)
          Times = np.array(Times)
          
     #class for pulling data files
     files = For_Plots(rundate)

     #Create lists of variable lists
     theta_file_list = [files.get_file(dump_time, "theta_bar") for dump_time in dump_time_list]
     press_file_list = [files.get_file(dump_time, "press") for dump_time in dump_time_list]
     flux_file_list = [files.get_file(dump_time, "wvelthetapert") for dump_time in dump_time_list]
     height_file = files.get_file("0000000600", "heights")

     AvProfLims = []
     invrinos = []
     gm_vars=[]
     #loop over text files files
     for i in range(len(theta_file_list)):
         
         theta = np.genfromtxt(theta_file_list[i])
         
         height = np.genfromtxt(height_file)
         t = Times[i]
         press = np.genfromtxt(press_file_list[i])
         rhow = nc.calc_rhow(press, height, theta[0])
         wvelthetapert = np.genfromtxt(flux_file_list[i])
         
         #only need up to 1900meters
         if rundate == "Jan152014_1":
             top_index = np.where(abs(2000 - height) < 26.)[0][0] #may need to be higher (e.g. for 60/2.5)
         else:
             top_index = np.where(abs(1700 - height) < 26.)[0][0] #may need to be higher (e.g. for 60/2.5)
             
         #function for calcuating heights
         [elbot_dthetadz, h, eltop_dthetadz, elbot_flux ,h_flux  ,eltop_flux, deltatheta, mltheta, z1_GM]= nc.Get_CBLHeights(height, press, theta, wvelthetapert, gamma, flux_s, top_index, 'old')
         [L0,N,B0,zenc]=nc.gm_vars(t,flux_s,gamma)
         
         
         delta_h=eltop_dthetadz - elbot_dthetadz
         delta = z1_GM - h
         
         print('z1_GM, delta: ',z1_GM, delta)
         
         [c_delta, rino, invrino, wstar, S, pi3, pi4] =  nc.calc_rino(B0, N, delta, h, zenc, mltheta, 1.0*flux_s/(rhow[0]*1004), deltatheta, gamma, delta_h)
         
         print(zenc, h, z1_GM, eltop_dthetadz, c_delta)
         
         AvProfLims.append([elbot_dthetadz, h, eltop_dthetadz, elbot_flux, h_flux, eltop_flux, deltatheta, mltheta, z1_GM])
         tau = 1.0*h/wstar
         invrinos.append([c_delta, rino, invrino, wstar, S, tau, mltheta, deltatheta, pi3, pi4])
         gm_vars.append([L0,N,B0,zenc])
     pdb.set_trace()
     files.save_file(np.array(AvProfLims), "AvProfLims_old")
     files.save_file(np.array(invrinos), "invrinos_old")
     #files.save_file(np.array(gm_vars), "gm_vars")

run_list = [ ["Mar12014", .01, 60],["Nov302013", .005, 100], ["Dec142013", .01, 100], ["Dec202013", .005, 60], ["Dec252013", .0025, 60], ["Jan152014_1", .005, 150], ["Mar52014", .01, 150]]

for run in run_list:
    Main_Fun(run[0], run[1], run[2])


    
    
