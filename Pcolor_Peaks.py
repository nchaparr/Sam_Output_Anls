from __future__ import division
from netCDF4 import Dataset
import glob,os.path
import numpy as np
import numpy.ma as ma
from scipy.interpolate import UnivariateSpline
from matplotlib import cm
from matplotlib import ticker
import matplotlib.pyplot as plt
#import site
#site.addsitedir('/tera/phil/nchaparr/SAM2/sam_main/python')
#from Percentiles import *
from matplotlib.patches import Patch
import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from Make_Timelist import *
import warnings
warnings.simplefilter('ignore', np.RankWarning)
#import pywt
from scipy import stats
from datetime import datetime

"""
    In testing phase -- get_fit() for identifying ML top
    To plot gradient maxima ie BL heights, and w on a 2d horizontal domain,
    and get a histogram or contour plot of BL heigths
    for an individual case
    added function to get ticks and labels based on mean and standard deviation
              
"""
#TODO: a mess right now.  but can be tidied up once regression code is included

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


def get_fit(theta, height):
     """
        Fitting the local theta profile with three lines
     
     """
     
     fitvals = np.zeros_like(theta)
     RSS = np.empty((298, 290))+ np.nan
     print RSS[0,0]
     for j in range(298):
          if j > 2:
               for k in range(298):
                    if k>j+1 and k<279:
                         b_1, a_1, num_b_11, num_b_12, num_b_13, dem_b_11, dem_b_12 = 0, 0, 0, 0, 0, 0, 0  
                         for i in range(j):
                              num_b_11 = num_b_11 + height[i]*theta[i]
                              num_b_12 = num_b_12 + height[i]
                              num_b_13 = num_b_13 + theta[i]
                              dem_b_11 = dem_b_11 + height[i]**2
                              dem_b_12 = dem_b_12 + height[i]
                         num_b = num_b_11 - 1/j*num_b_12*num_b_13
                         dem_b = dem_b_11  - 1/j*dem_b_12**2
                         b_1 = num_b/dem_b
                         a_1 = num_b_11/num_b_12 - b_1*dem_b_11/num_b_12

                         b_2, a_2, num_b_21, num_b_22, dem_b_21, dem_b_22, num_a_21, num_a_22 = 0, 0, 0, 0, 0, 0, 0, 0
                         for i in range(k-j):
                              num_b_21 = num_b_21 + theta[j+i]
                              num_b_22 = num_b_22 + a_1+b_1*height[j]
                              dem_b_21 = dem_b_21 + height[j+i]
                              dem_b_22 = dem_b_22 + height[j]
                              num_a_21 = num_a_21 + height[j]*theta[j]
                              num_a_22 = num_a_22 + height[j]**2
                         num_b2 = num_b_21 - num_b_22
                         dem_b2 = dem_b_21 - dem_b_22
                         b_2 = num_b2/dem_b2
                         a_2 = num_a_21/dem_b_22 - b_2*num_a_22/dem_b_21
                         
                         b_22 = (np.sum(theta[j:k]) - (k-j)*(a_1+b_1*height[j]))/(np.sum(height[j:k]) - (k-j)*height[j])                         
                         a_22 = np.sum(np.multiply(height[j:k], theta[j:k]))/np.sum(height[j:k]) - b_2*np.sum(height[j:k]**2)/np.sum(height[j:k])
                         
                         #print np.sum(theta[j:k]), num_b_21  NOTE: these are different, np.sum seems to round to 1 decimal place
                         
                         num_b_31, num_b_32, dem_b_31, dem_b_32, num_a_31, num_a_32 = 0, 0, 0, 0, 0, 0
                         for i in range(298-k):
                              num_b_31 = num_b_31 + theta[k+i]
                              num_b_32 = num_b_32 + a_2 + b_2*height[k]
                              dem_b_31 = dem_b_31 + height[k+i]
                              dem_b_32 = dem_b_32 + height[k]
                              num_a_31 = num_a_31 + height[k+i]*theta[k+i]
                              num_a_32 = num_a_32 + height[k+i]**2

                         num_b_3 = num_b_31 - num_b_32
                         dem_b_3 = dem_b_31 - dem_b_32

                         b_3 = num_b_3/dem_b_3
                         a_3 = num_a_31/dem_b_31 - b_2*num_a_32/dem_b_31
                                                  
                         #print np.sum(theta[k:298]), num_b_31 NOTE: sum seems to round up to 1 dec                                   
                         #print np.sum(np.multiply(height[k:298], theta[k:298])), num_a_31 NOTE: again rounding but this time not decimal places!!      
                                                                      
                         RSS[j, k] = np.sum(np.add(theta[2:j], -(a_1+ b_1*height[2:j]))**2) + np.sum(np.add(theta[j:k], -(a_2+ b_2*height[j:k]))**2) + np.sum(np.add(theta[k:298], -(a_2+ b_2*height[k:298]))**2) 
                         RSS_1 = 0
                         for i in range(j):
                              RSS_1 = RSS_1 + (theta[i] -(a_1 + b_1*height[i]))**2
                         RSS_2 = 0     
                         for i in range(k-j):
                              RSS_2 = RSS_2 + (theta[j+i] - (a_2 + b_2*height[j+i]))**2
                         RSS_3 = 0
                         for i in range(298-k):
                              RSS_3 = RSS_3 + (theta[k+i] - (a_3 + b_3*height[k+i]))**2
                         RSS[j, k] = RSS_1 + RSS_2 + RSS_3

                        
     RSS = ma.masked_where(np.isnan(RSS), RSS)
     [j, k] = np.unravel_index(ma.argmin(RSS), RSS.shape)
     
     b_1 = (np.sum(np.multiply(height[:j], theta[:j])) - 1/j*np.sum(height[:j]*np.sum(theta[:j])))/(np.sum(height[:j]**2) - 1/j*np.sum(height[2:j])**2)
     a_1 = np.sum(np.multiply(height[:j], theta[:j]))/np.sum(height[:j]) - b_1*np.sum(height[:j]**2)/np.sum(height[:j])
                         
     b_2 = (np.sum(theta[j:k]) - (k-j)*(a_1+b_1*height[j]))/(np.sum(height[j:k]) - (k-j)*height[j])                         
     a_2 = np.sum(np.multiply(height[j:k], theta[j:k]))/np.sum(height[j:k]) - b_2*np.sum(height[j:k]**2)/np.sum(height[j:k])

     b_3 = (np.sum(theta[k:298]) - (298-k)*(a_2+b_2*height[k]))/(np.sum(height[k:298]) - (298-k)*height[k])
     a_3 = np.sum(np.multiply(height[k:298], theta[k:298]))/np.sum(height[k:298]) - b_2*np.sum(height[k:298]**2)/np.sum(height[k:298])
                         
     
     fitvals[:j] = b_1*height[:j] + a_1
     fitvals[j:k] = b_2*height[j:k] + a_2
     fitvals[k:298] = b_3*height[k:298] + a_3


     return fitvals, RSS, j, k                                                              
                                                                                
                                                                              

#Lists of times relating to output (nc) files
dump_time_list, time_hrs = Make_Timelists(1, 600, 28800)
dump_time = dump_time_list[29]

for k in range(1):
     #getting variables from nc files
     [wvels, theta, tracer, height] = nc.Get_Var_Arrays("/tera2/nchaparr/Dec252013/runs/sam_case", "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_", dump_time, k+1)

     #getting points of maximum theta gradient, getting rid of this soon
     [dvardz, grad_peaks] = nc.Domain_Grad(theta, height) 
     tops_indices=np.where(np.abs(grad_peaks - 1400)<10)
     
     #choosing one horizontal point
     for i in range(4):
          top_index = [tops_indices[0][i], tops_indices[1][i]]
          [l, m] = top_index
          
          thetavals = theta[:, l, m]
          
          startTime = datetime.now()
          #print 'Start', startTime#1     
          fitvals, RSS, J, K = get_fit(thetavals, height)
          #print 'RSS time', (datetime.now()-startTime)
          #set up plot
          theFig = plt.figure(i)
          theFig.clf()
          theAx = theFig.add_subplot(121)
          theAx.set_title('')
          theAx.set_xlabel('')
          theAx.set_ylabel('z (m)')

          theAx1 = theFig.add_subplot(122)
          theAx1.set_title('')
          theAx1.set_xlabel('')
          theAx1.set_ylabel('')

          theAx1.plot(thetavals, height[:], 'wo')
          theAx.plot(fitvals[:J], height[:J], 'r-')
          theAx.plot(fitvals[J:K], height[J:K], 'b-')
          theAx.plot(fitvals[K:298], height[K:298], 'g-')
          print fitvals[:298].shape, height[:298].shape
          theAx1.plot(fitvals[:298], height[:298], 'r-')
          
          theAx1.set_xlim(300, 320)
          theAx1.set_ylim(0, 6000)
          theAx.set_ylim(0, 6000)
          theAx.set_xlim(300, 320)
          plt.show()




    
    
