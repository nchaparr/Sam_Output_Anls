from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import scipy.stats as stat
import matplotlib
import matplotlib.pyplot as plt
#import site
#site.addsitedir('/tera/phil/nchaparr/python/')
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from Make_Timelist import *

"""
  Gets 2d horizontal domain of gradient maxima ie BL heights
  
  Gets ELimits based on percentiles

  Calculates skewness of distributions.
  
"""

def Get_Var_Arrays(dump_time):
     """Pulls output from an ensemble cases, gets ensemble averages and perturbations and
     their horizontal averageses

    Arguments:
    dump_time -- time of output eg '0000000720'

    Returns:
    var_bar -- 64 array of horizontally averaged, ensemble averages or perturbations (covariances)
    
    """
     #create list of filenames for given dump_time
     ncfile_list = ["/tera2/nchaparr/Dec252013/runs/sam_case" + str(i+1) + "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_" + dump_time + ".nc" for i in range(9)]
     
     #create lists for variable arrays from each case
     wvels_list = []
     press_list = []
     height_list = []
     temp_list = []
     thetas_list = []
     tracer_list = []
     
     grad_max = []
     #ELLims_hist = []
     #ELLims_his1 = []
     for i in range(len(ncfile_list)): #loop over list of nc files
          thefile = ncfile_list[i]
          print thefile
          
          ncdata = Dataset(thefile,'r')
          wvel = np.squeeze(ncdata.variables['W'][...])
          
          press = np.squeeze(ncdata.variables['p'][...])#pressure already horizontally averaged
          height = np.genfromtxt('/tera/phil/nchaparr/python/Plotting/Dec252013/data/heights0000000600')
          temp = np.squeeze(ncdata.variables['TABS'][...])
          #tracer = np.squeeze(ncdata.variables['TRACER'][...])
          ncdata.close()
     
          #calculate thetas
          theta = np.zeros_like(temp)
          thetafact = np.array([(1.0*1000/k)**(1.0*287/1004) for k in press])
          for j in range(312):
               theta[j, :, :] = temp[j, :, :]*thetafact[j]

          thetas_list.append(theta)    #append array lists 
          wvels_list.append(wvel)
          
          press_list.append(press)
          height_list.append(height)
          temp_list.append(temp)
          #tracer_list.append(tracer)
          grad, peaks = nc.Domain_Grad(theta, height)
          grad_max.append(peaks)
     grad_max = np.array(grad_max)
     print grad_max.shape
     peaks = np.reshape(grad_max, (9*128*192,))     
     [ELBot, h, ELTop, Skew] = [np.percentile(peaks, 5), np.mean(peaks), np.percentile(peaks, 95), stat.skew(peaks, axis=0, bias=True)]
     #get arrays of enseble averaged variables
     ens_avthetas = nc.Ensemble1_Average(thetas_list)
     
     ens_avwvels = nc.Ensemble1_Average(wvels_list)
     ens_heights = height_list[0]
     ens_press = nc.Ensemble1_Average(press_list)
     #ens_tracer = nc.Ensemble1_Average(tracer_list)

     #for each point on the domain, get the gradients
     grad_ens_thetas, theta_peaks = nc.Domain_Grad(ens_avthetas, ens_heights)
     #find the point of maximum gradient above a certain level
     #ie clear of the surface
     #get max an min value for ELims and horizontal average is h
     theta_peaks = np.reshape(theta_peaks, 1*128*192)
     
     
     [ELBot1, h1, ELTop1, Skew1] = [np.min(theta_peaks), np.mean(theta_peaks), np.max(theta_peaks), stat.skew(theta_peaks, axis=0, bias=True)]   
     #append 
           
     return [ELBot, h, ELTop, Skew], [ELBot1, h1, ELTop1, Skew1] 

dump_time_list, time_hrs = Make_Timelists(1, 600, 28800)

if __name__ == "__main__":

     ELLims_hist = []
     ELLims_hist1 = []
     for dump_time in dump_time_list:
          [ELBot, h, ELTop, Skew], [ELBot1, h1, ELTop1, Skew1] = Get_Var_Arrays(dump_time)
          ELLims_hist.append([ELBot, h, ELTop, Skew])
          ELLims_hist1.append([ELBot1, h1, ELTop1, Skew1])
                              
     np.savetxt('/tera/phil/nchaparr/python/Plotting/Dec252013/data/ELLims_hist', np.array( ELLims_hist), delimiter=' ')
     np.savetxt('/tera/phil/nchaparr/python/Plotting/Dec252013/data/ELLims_hist1', np.array(ELLims_hist1), delimiter=' ')
     
          


    
    
