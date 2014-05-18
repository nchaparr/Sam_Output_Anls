from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
import matplotlib.pyplot as plt
#import site
#site.addsitedir('/tera/phil/nchaparr/python')
import nchap_fun as nc
from nchap_class import *
"""
   For contour plotting statistics output from an ensemble of runs
   Editing it for plotting scaled TKE   
"""
#TODO: may be obsolete
def Get_Var_Arrays(var, fignum):
     """Pulls stats output from a ensemble of cases, gets ensemble averages and does contour plots
     on height, time grid

    Arguments:
    var -- key in ''
    fignum  -- integter for figure
     
    
    """
     #create list of filenames
     ncfile_list = ["/tera2/nchaparr/Mar52014/runs/sam_case"+ str(i+1) + "/OUT_STAT/NCHAPP1_testing_doscamiopdata.nc" for i in range(10)]

     #create lists for variable arrays from each case
     vars_list = []
     height_list = []
     press_list = []
     time_list = []
     
     for i in range(len(ncfile_list)): #loop over list of nc files
          thefile = ncfile_list[i]
          ncdata = Dataset(thefile,'r')
          Vars = ncdata.variables[var][...]
          #print Vars.shape
          press = ncdata.variables['PRES'][...]
          height = ncdata.variables['z'][...]
          top = np.where(abs(height - 2000) < 50)[0][0]
          Vars = Vars[:, 0:top]
          height = height[0:top]
          time = ncdata.variables['time'][...]
          ncdata.close()
               
          vars_list.append(Vars)          
          height_list.append(height)
          time_list.append(time)    
          press_list.append(press)
          
     #get ensemble averages
     ens_vars = nc.Ensemble1_Average(vars_list)
     ens_press = nc.Ensemble1_Average(press_list)
     ens_press = np.transpose(ens_press)
     
     #TODO: verify this is in time order!
     
     #print 'ENSEMBLE AVERAGED',  ens_vars.shape
     time = (time_list[0]-169)*24
     height = height_list[0] #TODO: time, height don't need to be averaged 
          
     #set up plot
     theAx = nc.Do_Plot(fignum, 'Scaled ' + var + ' vs Scaled Height', 'Height/h', var+'/w*2', 111)
     #print ens_vars.shape, height.shape
     have_ens_vars = [] 
     print len(time)
     for i in range(len(time)):
          if np.mod(i+1, 6)==0:
               #print i, time[i], 1.0*(i+1)/10, "plotting"
               points = For_Plots("Mar52014")
               rinovals = points.rinovals()
               print len(rinovals[:,2])
               AvProfVars = points.AvProfVars()
               #invrinos: [rino, invrino, wstar, S, tau, mltheta, deltatheta, pi3, pi4]
               wstar= rinovals[1.0*((i+1))*(6.0/6.0)-1, 2]
               h= AvProfVars[1.0*((i+1))*(6.0/6.0)-1, 1]
               h_index = np.where(height==h)[0]
               print time[i]
               have_ens_vars.append(1.0*np.sum(ens_vars[i][0:h_index])/(h*wstar**3))
               #print have_ens_vars[i]
               theAx.plot(1.0*ens_vars[i]/wstar**3, 1.0*height/h, label=str(int(time[i])+1) + 'hrs')
     #height, time = np.meshgrid(height, time)
     #maxlev = np.max(ens_vars)
     #minlev = np.min(ens_vars)
     #step = (maxlev- minlev)/20
     #levels = [i for i in np.arange(minlev, maxlev, step)]
     #CS = plt.contourf(time, height, ens_vars, levels, cmap=plt.cm.bone)
     #cbar = plt.colorbar(CS)
     #cbar.ax.set_ylabel('colorbar')
     #print 'plotting'
     #theAx.plot(time, have_ens_vars, label=var)
     plt.ylim(0, 2)
     plt.legend(loc = 'upper right', prop={'size':8})
     plt.show()

#theAx = nc.Do_Plot(1, 'Layer Averaged, Scaled TKE Terms vs Time', 'TKE Term/w*3', 'Time (hrs)',111)
var_list = [ 'TKE']
#BUOYA', 'BUOYAS', 'DISSIP', 'DISSIPS','DISSIPS', 'BUOYAS', TKE', 'TKES','TKE', 'TKES', 'WVADV', 'WUADV', 'WUPRES', 'WVPRES', 'WUSHEAR', 'WVSHEAR', 'W2ADV', 'W2PRES', 'W2BUOY', 'WVBUOY', 'WUBUOY', 'W2REDIS', 'W2DIFF'
for i in range(len(var_list)):
     Get_Var_Arrays(var_list[i], i)
#plt.ylim()
#plt.ylim(-.000035, .000035)
#plt.legend(loc='lower right')
#plt.show()    



    
    
