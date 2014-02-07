from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
import matplotlib.pyplot as plt
import site
#site.addsitedir('/tera/phil/nchaparr/SAM2/sam_main/python')
#from nchap_fun import *
from Make_Timelist import *
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc


"""
   For comparing 3D output to statistics output
      
"""

def Get_Var_Arrays(dump_time):
     """Pulls output nstcdf file

    Arguments:
    dump_time -- time of output eg '0000000720'

    Returns:
    var_bar -- 64 array of horizontally averaged, ensemble averages or perturbations (covariances)
    
    """

      #create filename for given dump_time
     ncfile = "/tera/phil/nchaparr/sam_ensemble/sam_case"+ str(1) + "/OUT_3D/TOGA_1_testing_doscamiopdata_16_" + dump_time + ".nc"
     
     #loop over list of nc files
     thefile = ncfile
     print thefile
     ncdata = Dataset(thefile,'r')
     wvel = np.squeeze(ncdata.variables['W'][...])
     press = np.squeeze(ncdata.variables['p'][...])#pressure already horizontally averaged
     height = np.squeeze(ncdata.variables['z'][...])
     temp = np.squeeze(ncdata.variables['TABS'][...])
     
     ncdata.close()
     
          #calculate thetas
     theta = np.zeros([64, 256, 256])
     thetafact = np.array([(1.0*1000/k)**(1.0*287/1004) for k in press])
     for j in range(64):
          theta[j, :, :] = temp[j, :, :]*thetafact[j]

    #get horizontally averaged wprimethetaprimes and plot
     theta_bar = nc.Horizontal_Average(theta)

     rows, cols, cols1 = theta.shape

     #print rows, cols, cols1, theta.shape

     theta_pert = np.zeros_like(theta)

     for row in range(rows):
          if row == 0:
               theta_pert[row] = (theta[row, :, :] - theta_bar[row])
          else:
               theta_pert[row] = 0.5*(theta[row, :, :] - theta_bar[row] + theta[row-1, :, :] - theta_bar[row-1])
               
     wveltheta = np.multiply(theta_pert, wvel)
     wveltheta_bar = nc.Horizontal_Average(wveltheta)
     print'shape of theta_pert', theta_pert.shape
     theta_pert_bar = nc.Horizontal_Average(theta_pert)     
     #print 'height_bar', height_bar

     return wveltheta_bar, height

dump_time_list, Times = Make_Timelists(2, 30, 14400)

#set up plot
theFig = plt.figure(1)
theFig.clf()
theAx = theFig.add_subplot(111)
theAx.set_title('')
theAx.set_xlabel('')
theAx.set_ylabel('Height (m)')
#get horizontally averaged wprimethetaprimes and plot
for i in range(len(dump_time_list)):     
     #hor_avs = Horizontal_Average(wvelthetaperts_list[i])
     #print len(hor_avs), len(height_bar), hor_avs, height_bar
     wvelpert_bar, height_bar = Get_Var_Arrays(dump_time_list[i])
     wvelpert_bar[0] = np.nan
     #wvelpert_bar = np.multiply(wvelpert_bar, np.zeros_like(wvelpert_bar)+1004)
     theAx.plot(wvelpert_bar, height_bar, label = dump_time_list[i])
     zeros = np.zeros(len(height_bar))
theAx.plot(zeros, height_bar)
plt.legend(loc = 'upper right', prop={'size':8})
plt.ylim(0, 2500)
#plt.xlim(295, 310)
plt.show()



    
    
