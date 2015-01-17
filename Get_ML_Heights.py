from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker
import site
site.addsitedir('/tera2/nchaparr/Dec252013/python')
import sys
#sys.path.insert(0, '/newtera/teraphil/nchaparr/python')
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
def Main_Fun(dump_time_index, case, date):
     """Loops over ensemble cases.  Pulls temperature, pressure, height from nc output files using nchap_class
    gets height of mixed layer at each horizontal point using fast_fit and saves as a txt file.   

    Arguments:
    dump_time -- time of output eg '0000000720'
    case -- integer

    Returns:
    ML_Heights -- 128*192 array as txt file
    
    """
     dump_time_list, time_hrs = Make_Timelists(1, 3600, 28800)
    
     ncfile = "/newtera/tera/phil/nchaparr/tera2_cp/nchaparr/"+ date +"/runs/sam_case" + str(case) + "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_" + dump_time_list[dump_time_index] + ".nc"
               
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
     for i in range(128): #change to 128          
          for j in range(192):#change to 192
               top = np.where(np.abs(height-2300)<100)[0][0]
               RSS, J, K = fsft.get_fit(theta[:, i, j], height, top)
               print i, j, J
               ML_Heights[i, j] = height[J]
               #if height[J] < 500:
               #print i, j, J, height[J]#TODO:this could definitely be a function
               b_1 = (np.sum(np.multiply(height[:J], theta[:J, i, j])) - 1.0/J*np.sum(height[:J])*np.sum(theta[:J, i, j]))/(np.sum(height[:J]**2) - 1.0/J*np.sum(height[:J])**2)
               a_1 = np.sum(np.multiply(height[:J], theta[:J, i, j]))/np.sum(height[:J]) - b_1*np.sum(height[:J]**2)/np.sum(height[:J])                          
                                   
               b_2 = (np.sum(theta[J:K, i, j]) - (K-J)*(a_1+b_1*height[J]))/(np.sum(height[J:K]) - (K-J)*height[J])                                    
               a_2 = np.sum(np.multiply(height[J:K], theta[J:K, i, j]))/np.sum(height[J:K]) - b_2*np.sum(height[J:K]**2)/np.sum(height[J:K])
               b_3 = (np.sum(theta[K:290, i, j]) - (290-K)*(a_2+b_2*height[K]))/(np.sum(height[K:290]) - (290-K)*height[K])               
               if b_3-b_2>.0017:
                   #print i, j, J, height[J], height[K], b_3, b_2, b_3-b_2
                    ML_Heights[i, j] = height[K]
                    
     print '/newtera/tera/phil/nchaparr/python/Plotting/'+date+'/data/mixed_layer_height_'+ str(case) + '_' + dump_time_list[dump_time_index]           
     np.savetxt('/newtera/tera/phil/nchaparr/python/Plotting/'+date+'/data/mixed_layer_height_'+ str(case) + '_' + dump_time_list[dump_time_index], ML_Heights, delimiter=' ')
     
     return ML_Heights

def Call_Main_Fun(date):
     dump_time_no = 8
     case_no = 10
     date = 'Dec252013'
     for i in range(dump_time_no):
          if i == 4:
               for j in range(case_no):
                    ML_Heights = Main_Fun(i, j+1, date)
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
         
                    im = theAx.pcolor(X, Y, np.transpose(ML_Heights), cmap=cm.hot, vmax=v_max, vmin=v_min)
                    bar = plt.colorbar(im)
                    plt.xlim(0, 3200)
                    plt.ylim(0, 4800)
               
                    #label_list, tick_list = get_ticks(mean, stddev,v_max, v_min)
     
                    #bar.locator = ticker.FixedLocator(tick_list)
                    #bar.formatter= ticker.FixedFormatter(label_list)
                    #bar.update_ticks()
                    plt.show()
                    #print time_hrs[i]
                    #theFig.savefig('/newtera/tera/phil/nchaparr/python/Plotting/'+date+'/pngs/h_cont_'+str(time_hrs[i])+'_hrs.png')
     

dump_time_list, time_hrs = Make_Timelists(1, 3600, 28800)
print time_hrs
date = "Dec252013"

if __name__ == "__main__":

     #import argparse
     #parser = argparse.ArgumentParser()
     #parser.add_argument("--run_name", type=str, help='name of run')
     #args = parser.parse_args()
     #date = args.run_name
     #print date
     Call_Main_Fun(date)
     
     #set up plot



    
    
