from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
import matplotlib.pyplot as plt
#import site
#site.addsitedir('/tera/phil/nchaparr/python/')
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from Make_Timelist import *
import fastfit as fsft

"""
  Gets 2d horizontal domain of gradient maxima ie BL heights
  and plots histograms thereof.
      
"""
#TODO: needs to change not that maximum gradient isn't to be used
def Main_Fun(dump_time):
     """Loops over ensemble cases.  Pulls temperature, pressure, height from nc output files using nchap_class
    gets height of mixed layer at each horizontal point using fast_fit and saves as a txt file.   

    Arguments:
    dump_time -- time of output eg '0000000720'

    Returns:
    ML_Heights -- 
    
    """
     #create list of filenames for given dump_time
     ncfile_list = ["/tera2/nchaparr/Dec252013/runs/sam_case" + str(i+1) + "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_" + dump_time + ".nc" for i in range(9)]
     #print ncfile_list
     #create lists for variable arrays from each case  
     
     for i in range(len(ncfile_list)): #loop over list of nc files
          thefile = ncfile_list[i]
          print thefile
          ncdata = Dataset(thefile,'r')
          #wvel = np.squeeze(ncdata.variables['W'][...])
          
          press = np.squeeze(ncdata.variables['p'][...])#pressure already horizontally averaged
          height = np.squeeze(ncdata.variables['z'][...])
          temp = np.squeeze(ncdata.variables['TABS'][...])
          #tracer = np.squeeze(ncdata.variables['TRACER'][...])
          ncdata.close()
     
          #calculate thetas
          theta = np.zeros_like(temp)
          thetafact = np.array([(1.0*1000/k)**(1.0*287/1004) for k in press])
          for j in range(312):
               theta[j, :, :] = temp[j, :, :]*thetafact[j]

          ML_Heights = np.empty([128, 192])
          for i in range(128):
               for j in range(192):
                    RSS, J, K = fsft.fastfit(theta[:, i, j], height)
                    ML_Heights[i, j] = height[J]          
          np.savetxt('/tera/phil/nchaparr/python/Plotting/Dec252013/data/mixed_layer_height_'+ str(i+1) + '_' + dump_time, ML_Heights, delimiter=' ')



     peaks = np.reshape(grad_max, (1, 9*128*192))     
     #get arrays of enseble averaged variables
     ens_avthetas = nc.Ensemble1_Average(thetas_list)
     #print 'ens ave thetas', ens_avthetas
     ens_avwvels = nc.Ensemble1_Average(wvels_list)
     ens_heights = height_list[0]
     ens_press = nc.Ensemble1_Average(press_list)
     #ens_tracer = nc.Ensemble1_Average(tracer_list)

     #for each point on the domain, get the gradients
     ens_grad, ens_peaks = nc.Domain_Grad(ens_avthetas, ens_heights)
     #find the point of maximum gradient above a certain level
     #ie clear of the surface
     
     #do a histogram         
     
     return ens_peaks, ens_heights, peaks

dump_time_list, time_hrs = Make_Timelists(1, 600, 28800)

if __name__ == "__main__":
 
     #set up plot
     theFig = plt.figure(1)
     theFig.clf()
     theAx = theFig.add_subplot(111)
     theAx.set_title(r"$Histogram \ of \ Local \ h \ after \ 1 \ hours$")
     theAx.set_xlabel(r"$z \ (m)$")
     theAx.set_ylabel(r"$Number \ in \ z \ Bin$")
     dump_time = dump_time_list[5]
     [ens_peaks, heights, peaks] = Main_Fun(dump_time)
     #n, bins, patches = theAx.hist(peaks, bins=20)
     height_bin_vols = nc.Bin_Peaks(peaks, heights)
     theAx.bar(heights, height_bin_vols)
     #plt.legend(loc = 'upper left', prop={'size':8})
     plt.ylim(0, 25000)
     plt.xlim(0, 2000)
     plt.show()
     



    
    
