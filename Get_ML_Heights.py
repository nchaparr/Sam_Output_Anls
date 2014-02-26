from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker
#import site
#site.addsitedir('/tera/phil/nchaparr/python/')
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from Make_Timelist import *
import fastfit as fsft

"""
  Gets 2d horizontal domain of ML heights
  and plots pcolor plots thereof.
      
"""

def get_ticks(mean, stddev, max, min):

     """

     gets ticks and tick lavels for contour plot based on mean and standard deviation
    
     Arguments:    
     mean, stddev, max, min

     Returns:       
     ticks, tick_labels
    
     """
     tick_list = []
     label_list = []     
     int1=int(np.ceil((mean-min)/stddev))     
     int2=int(np.ceil((max-mean)/stddev))
     
     
     
     for i in range(int1):
         if int1==1:
             tick_list.append(min)
             label_list.append(r'$\mu - %.1f \sigma$' %((mean-min)/stddev))
             
         elif i > 0:
             tick_list.append(mean - (int1-i)*stddev)
             
             label_list.append(r'$\mu - %.1f \sigma$' %(int1-i))
             
         #else:
             #tick_list.append(min)
             
             #label_list.append(r'$\mu - %.1f \sigma$' %((mean-min)/stddev))
             
     tick_list.append(mean)       
     label_list.append(r'$\mu$')
     
     
     for i in range(int2):
         
         if int2==1:
             tick_list.append(max)
             
             label_list.append(r'$\mu + %.1f \sigma$' %((max-mean)/stddev))
             
         elif i< int2-1:
             tick_list.append(mean + (i+1)*stddev)
             
             label_list.append(r'$\mu + %.1f \sigma$' %(i+1))
             
         #else:
             #tick_list.append(max)
             
             #label_list.append(r'$\mu + %.1f \sigma$' %((max-mean)/stddev))
             
     return label_list, tick_list


#TODO: needs to change not that maximum gradient isn't to be used
def Main_Fun(dump_time, case):
     """Loops over ensemble cases.  Pulls temperature, pressure, height from nc output files using nchap_class
    gets height of mixed layer at each horizontal point using fast_fit and saves as a txt file.   

    Arguments:
    dump_time -- time of output eg '0000000720'
    case -- integer

    Returns:
    ML_Heights -- 128*192 array as txt file
    
    """
     #create list of filenames for given dump_time
     ncfile = "/tera2/nchaparr/Dec202013/runs/sam_case" + str(case) + "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_" + dump_time + ".nc"
     
          
     thefile = ncfile
     print thefile
     ncdata = Dataset(thefile,'r')
          
          
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
               RSS, J, K = fsft.get_fit(theta[:, i, j], height)
               ML_Heights[i, j] = height[J]
               if height[J] < 500:
                    print i, j, J, height[J]
               
     print '/tera/phil/nchaparr/python/Plotting/Dec202013/data/mixed_layer_height_'+ str(case) + '_' + dump_time           
     np.savetxt('/tera/phil/nchaparr/python/Plotting/Dec202013/data/mixed_layer_height_'+ str(case) + '_' + dump_time, ML_Heights, delimiter=' ')
     
     return ML_Heights

dump_time_list, time_hrs = Make_Timelists(1, 3600, 28800)

if __name__ == "__main__":
 
     #set up plot
     
     for i in range(len(dump_time_list)):
          if i == 4:
               for j in range(1):
                    ML_Heights = Main_Fun(dump_time_list[i], j+1)

               theFig = plt.figure(i)
               theFig.clf()
               theAx = theFig.add_subplot(111)
               theAx.set_title(r"$Contour \ of \ Local \ h \ after \ " + str(time_hrs[i]) +" \ hours$")
               theAx.set_xlabel(r"$x \ (m)$")
               theAx.set_ylabel(r"$y \ (m)$")

               v_max, v_min, mean, stddev = np.amax(ML_Heights), np.amin(ML_Heights), np.mean(ML_Heights), np.std(ML_Heights)

               filler_array = np.zeros([64, 192])
               ML_Heights = np.vstack((ML_Heights, filler_array))
               x = np.arange(0, 4800, 25)
               y = np.arange(0, 4800, 25)
               X,Y = np.meshgrid(x, y)
         
               im = theAx.pcolor(X, Y, ML_Heights, cmap=cm.bone, vmax=v_max, vmin=v_min)
               bar = plt.colorbar(im)
               plt.ylim(0, 3200)
               plt.xlim(0, 4800)
               
               #v_max, v_min, mean, stddev = np.amax(grad_peaks), np.amin(grad_peaks), np.mean(grad_peaks), np.std(grad_peaks)
               label_list, tick_list = get_ticks(mean, stddev,v_max, v_min)
     
               bar.locator = ticker.FixedLocator(tick_list)
               bar.formatter= ticker.FixedFormatter(label_list)
               bar.update_ticks()
               plt.show()
               print time_hrs[i]
               theFig.savefig('/tera/phil/nchaparr/python/Plotting/Dec202013/pngs/cont_'+str(time_hrs[i])+'_hrs.png')
     



    
    
