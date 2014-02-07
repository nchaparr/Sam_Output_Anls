from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
import matplotlib.pyplot as plt
import site
site.addsitedir('/tera/phil/nchaparr/python')
import nchap_fun as nc

"""
   For contour plotting statistics output from an ensemble of runs
      
"""
#TODO: may be obsolete
#TODO: see which functions can be pulled from nchap_fun
def Ensemble_Averageold(list, dim2):
     """Gets enseble average of a list of arrays

    Arguments:
    list -- list of 6, dim2 arrays

    Returns:
    ens_avs -- 6, dim2 array
    
    """
     ens_avs = np.zeros([24, dim2])
     for i in range(24): #time
          for j in range(dim2): #height
               vals = []          
               for k in range(len(list)):
                    #print i, j, k, list[k].shape
                    val = list[k][i, j]#i, j, kth elemement from var array l                                        
                    vals.append(val)                
               avval = 1.0*sum(vals)/len(vals) #average over l arrays                             
               ens_avs[i, j] = avval
     return ens_avs

def Ensemble_Average(list):
    """Gets enseble average of a list of arrays

    Arguments:
    list -- list of arrays

    Returns:
    ens_avs -- array
    
    """
    to_av = list[0]
    for k in range(len(list)-1): 
        to_av = np.add(to_av, list[k+1])
    ens_avs = 1.0*to_av/len(list)
    return ens_avs


def Get_Var_Arrays(var, fignum):
     """Pulls stats output from a ensemble of cases, gets ensemble averages and does contour plots
     on height, time grid

    Arguments:
    var -- key in ''
    fignum  -- integter for figure
     
    
    """
     #create list of filenames
     ncfile_list = ["/tera/phil/nchaparr/python/Plotting/Sep302013/case"+ str(i+1) + "/OUT_STAT/NCHAPP1_testing_doscamiopdataVI.nc" for i in range(10)]

     #create lists for variable arrays from each case
     vars_list = []
     height_list = []
     press_list = []
     time_list = []
     
     for i in range(len(ncfile_list)): #loop over list of nc files
          thefile = ncfile_list[i]
          ncdata = Dataset(thefile,'r')
          Vars = ncdata.variables[var][...]
          press = ncdata.variables['PRES'][...]
          height = ncdata.variables['z'][...]
          top = np.where(abs(height - 2500) < 50)[0][0]
          Vars = Vars[:, 0:top]
          height = height[0:top]
          time = ncdata.variables['time'][...]
          ncdata.close()
               
          vars_list.append(Vars)          
          height_list.append(height)
          time_list.append(time)    
          press_list.append(press)
          
     #get ensemble averages
     ens_vars = Ensemble_Average(vars_list)
     ens_press = Ensemble_Average(press_list)
     ens_press = np.transpose(ens_press)
     
     #TODO: verify this is in time order!
     
     #print 'ENSEMBLE AVERAGED',  ens_vars.shape
     time = time_list[0]
     height = height_list[0] #TODO: time, height don't need to be averaged 
          
     #set up plot
     theAx = nc.Do_Plot(fignum, var + ' vs time vs height', 'Height (m)', 'time', 111)
     height, time = np.meshgrid(height, time)
     maxlev = np.max(ens_vars)
     minlev = np.min(ens_vars)
     step = (maxlev- minlev)/20
     levels = [i for i in np.arange(minlev, maxlev, step)]
     CS = plt.contourf(time, height, ens_vars, levels, cmap=plt.cm.bone)
     cbar = plt.colorbar(CS)
     cbar.ax.set_ylabel('colorbar')
     plt.show()

var_list = ['THETAV', 'TVFLUX']
#'TKE', 'TKES', 'WVADV', 'WUADV', 'WUPRES', 'WVPRES', 'WUSHEAR', 'WVSHEAR', 'W2ADV', 'W2PRES', 'W2BUOY', 'WVBUOY', 'WUBUOY', 'W2REDIS', 'W2DIFF'
for i in range(len(var_list)):
     Get_Var_Arrays(var_list[i], i)




    
    
