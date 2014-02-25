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
import fastfit as fsft

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
     RSS = np.empty((290, 290))+ np.nan
     print RSS[0,0]
     for j in range(290):
          if j > 2:
               for k in range(290):
                    if k>j+1 and k<289:
                         b_1 = (np.sum(np.multiply(height[:j], theta[:j])) - 1/j*np.sum(height[:j])*np.sum(theta[:j]))/(np.sum(height[:j]**2) - 1/j*np.sum(height[:j])**2)
                         a_1 = np.sum(np.multiply(height[:j], theta[:j]))/np.sum(height[:j]) - b_1*np.sum(height[:j]**2)/np.sum(height[:j])
                         
                         b_2 = (np.sum(theta[j:k]) - (k-j)*(a_1+b_1*height[j]))/(np.sum(height[j:k]) - (k-j)*height[j])
                         
                         a_2 = np.sum(np.multiply(height[j:k], theta[j:k]))/np.sum(height[j:k]) - b_2*np.sum(height[j:k]**2)/np.sum(height[j:k])

                         b_3 = (np.sum(theta[k:290]) - (290-k)*(a_2+b_2*height[k]))/(np.sum(height[k:290]) - (290-k)*height[k])
                         a_3 = np.sum(np.multiply(height[k:290], theta[k:290]))/np.sum(height[k:290]) - b_3*np.sum(height[k:290]**2)/np.sum(height[k:290])
                         
                         RSS[j, k] = np.sum(np.add(theta[2:j], -(a_1+ b_1*height[2:j]))**2) + np.sum(np.add(theta[j:k], -(a_2+ b_2*height[j:k]))**2) + np.sum(np.add(theta[k:290], -(a_3+ b_3*height[k:290]))**2) 
                         
                         
     RSS = ma.masked_where(np.isnan(RSS), RSS)
     [j, k] = np.unravel_index(ma.argmin(RSS), RSS.shape)
     
     b_1 = (np.sum(np.multiply(height[:j], theta[:j])) - 1/j*np.sum(height[:j]*np.sum(theta[:j])))/(np.sum(height[:j]**2) - 1/j*np.sum(height[2:j])**2)
     a_1 = np.sum(np.multiply(height[:j], theta[:j]))/np.sum(height[:j]) - b_1*np.sum(height[:j]**2)/np.sum(height[:j])
                         
     b_2 = (np.sum(theta[j:k]) - (k-j)*(a_1+b_1*height[j]))/(np.sum(height[j:k]) - (k-j)*height[j])                         
     a_2 = np.sum(np.multiply(height[j:k], theta[j:k]))/np.sum(height[j:k]) - b_2*np.sum(height[j:k]**2)/np.sum(height[j:k])

     b_3 = (np.sum(theta[k:290]) - (290-k)*(a_2+b_2*height[k]))/(np.sum(height[k:290]) - (290-k)*height[k])
     a_3 = np.sum(np.multiply(height[k:290], theta[k:290]))/np.sum(height[k:290]) - b_3*np.sum(height[k:290]**2)/np.sum(height[k:290])
                         
     
     fitvals[:j] = b_1*height[:j] + a_1
     fitvals[j:k] = b_2*height[j:k] + a_2
     fitvals[k:290] = b_3*height[k:290] + a_3


     return fitvals, RSS, j, k                                                
                                                                                
                                                                              

#Lists of times relating to output (nc) files
dump_time_list, time_hrs = Make_Timelists(1, 600, 28800)
dump_time = dump_time_list[29]
print dump_time

for k in range(1):
     #getting variables from nc files
     [wvels, theta, tracer, height] = nc.Get_Var_Arrays("/tera2/nchaparr/Dec252013/runs/sam_case", "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_", dump_time, k+1)

     #getting points of maximum theta gradient, getting rid of this soon
     #[dvardz, grad_peaks] = nc.Domain_Grad(theta, height) 
     #tops_indices=np.where(np.abs(grad_peaks - 1400)<10)
     
     #choosing one horizontal point
     for i in range(1):
          #top_index = [tops_indices[0][i], tops_indices[1][i]]
          #[i, j] = top_index
          [i, j] = [53, 151]
          thetavals = theta[:, i, j]

          startTime = datetime.now()
          print 'Start', startTime#1     
          RSS, J, K = fsft.get_fit(thetavals, height)
          print J, height[J]
          print 'RSS time', (datetime.now()-startTime)
          fitvals = np.zeros_like(thetavals)
          b_1 = (np.sum(np.multiply(height[9:J], thetavals[9:J])) - 1.0/(J-9)*np.sum(height[9:J]*np.sum(thetavals[9:J])))/(np.sum(height[9:J]**2) - 1.0/(J-9)*np.sum(height[9:J])**2)
          print np.sum(np.multiply(height[9:J], thetavals[9:J])),  - 1.0/(J-9)*np.sum(height[9:J]*np.sum(thetavals[9:J])), np.sum(height[9:J]**2), - 1.0/(J-9)*np.sum(height[9:J])**2

          a_1 = np.sum(np.multiply(height[9:J], thetavals[9:J]))/np.sum(height[9:J]) - b_1*np.sum(height[9:J]**2)/np.sum(height[9:J])
                         
          b_2 = (np.sum(thetavals[J:K]) - (K-J)*(a_1+b_1*height[J]))/(np.sum(height[J:K]) - (K-J)*height[J])                         
          a_2 = np.sum(np.multiply(height[J:K], thetavals[J:K]))/np.sum(height[J:K]) - b_2*np.sum(height[J:K]**2)/np.sum(height[J:K])

          b_3 = (np.sum(thetavals[K:290]) - (290-K)*(a_2+b_2*height[K]))/(np.sum(height[K:290]) - (290-K)*height[K])
          a_3 = np.sum(np.multiply(height[K:290], thetavals[K:290]))/np.sum(height[K:290]) - b_3*np.sum(height[K:290]**2)/np.sum(height[K:290])
                         
     
          fitvals[:J] = b_1*height[:J] + a_1
          fitvals[J:K] = b_2*height[J:K] + a_2
          fitvals[K:290] = b_3*height[K:290] + a_3

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
          theAx.plot(fitvals[K:290], height[K:290], 'g-')
          theAx1.plot(fitvals[:290], height[:290], 'r-')

theAx1.set_xlim(300, 320)
theAx1.set_ylim(0, 4000)
theAx.set_ylim(0, 4000)
theAx.set_xlim(300, 320)
plt.show()




    
    
