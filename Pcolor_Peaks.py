from __future__ import division
from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib import cm
from matplotlib import ticker
import matplotlib.pyplot as plt
#import site
#site.addsitedir('/tera/phil/nchaparr/SAM2/sam_main/python')
#from Percentiles import *
from matplotlib.patches import Patch
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from Make_Timelist import *
import warnings
warnings.simplefilter('ignore', np.RankWarning)
#import pywt



"""
    now loops over all cases at one time
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


def get_fit(array, height):

    """
        This is where I'll put the regression code

    Arguments:
    --  array

    Returns:
     -- array    
    
    """
    from scipy.interpolate import splrep, splev, PiecewisePolynomial
    
    
    
    [znum, ynum, xnum] = array.shape
    chunk_list=[]
    count=0
      
    while ((count+1)*50 + count*10)<znum:         
         bottom = count*50 - count*10
         top = (count+1)*50
         print count, bottom, top
         chunk_list.append([bottom, top])
         count = count + 1
    chunk_list.append([top-10, znum])

    print count, top-10, znum

    fit_array = np.zeros_like(array)
    grad_array = np.zeros_like(array)
    for i in range(ynum):
         for j in range(xnum):
              for k in range(len(chunk_list)):
                   
                   bottom = chunk_list[k][0]
                   top = chunk_list[k][1]
                   
                   #fitfunc = np.polyfit(height[bottom:top], array[bottom:top, i, j], 3)
                #   if k == 0:
                        #print k
                        #fit = np.add(fitfunc[0]*height[bottom:top-10]**3, fitfunc[1]*height[bottom:top-10]**2)
                        #print 'fit array shape', fit.shape
                        #fit = np.add(fit, fitfunc[2]*height[bottom:top-10]**1) + fitfunc[3]
                        
                        #print 'fit_array shape', fit.shape
                        #fit = np.add(fit, fitfunc[3]*height[bottom:top-5]**1) + fitfunc[4]
                        
               #         grad = np.add(3*fitfunc[0]*height[bottom:top-10]**2, 2*fitfunc[1]*height[bottom:top-10]**1) + fitfunc[2] 
                        #grad = np.add(grad, 2*fitfunc[2]*height[bottom:top-5]) + fitfunc[3]
                        #print grad.shape, fit.shape
                        #print 'array shapes', fit.shape, grad.shape
                   
              #     else:
                        #print k
                        #print 'array shapes', fit.shape, grad.shape
                       # fit1 = np.add(fitfunc[0]*height[bottom+10:top-10]**3, fitfunc[1]*height[bottom+10:top-10]**2)
                       # fit1 = np.add(fit1, fitfunc[2]*height[bottom+10:top-10]**1) + fitfunc[3]
                       # grad1 = np.add(3*fitfunc[0]*height[bottom+10:top-10]**2, 2*fitfunc[1]*height[bottom+10:top-10]**1) + fitfunc[2] 
                       # fit = np.hstack((fit, fit1))
                       # grad = np.hstack((grad, grad1))
                       # print grad.shape, fit.shape

              tck = splrep(height[10:], array[10:, i, j], s=0)
              #print tck
              fit = splev(height[10:], tck, der=0)
              fit = np.hstack((np.zeros([10,])+fit[0], fit))
              grad = splev(height[10:], tck, der=2)
              grad = np.hstack((np.zeros([10,]), grad))
              fit_array[:, i, j] = fit              
              grad_array[:, i, j] = grad              
              
    grad_peaks = np.zeros((ynum, xnum))
   
    
    bot_index = np.where(abs(250 - height) < 5)[0][0]
    top_index = np.where(abs(1700 - height) < 5.)[0][0]
    for i in range(ynum):
         for j in range(xnum):               
              max_grad = np.amax(grad_array[bot_index:top_index, i, j])               
              index = np.where(grad_array[bot_index:top_index, i, j] - max_grad == 0)[0][0]
              print height[bot_index], height[top_index], height[index+bot_index]
              BLheight = height[index+bot_index]               
              grad_peaks[i, j] = BLheight          
     
    return fit_array, grad_array, grad_peaks


def get_RSS(theta, height):
     """
        Fitting the local theta profile with three lines
     
     """
     fitvals = np.zeros_like(theta)
     RSS = np.empty((312, 312))+ np.nan
     for j in range(312):
          if j > 2:
               for k in range(312):
                    if k>j+1:
                         b_1 = (np.sum(np.multiply(height[2:j], theta[2:j])) - 1/j*np.sum(height[2:j]*np.sum(theta[2:j])))/(np.sum(height[2:j]**2) - 1/j*np.sum(height[2:j])**2)
                         a_1 = np.sum(np.multiply(height[2:j], theta[2:j]))/np.sum(height[2:j]) - b_1*np.sum(height[2:j]**2)/np.sum(height[2:j])
                         
                         b_2 = (np.sum(theta[j:k]) - (k-j)*(a_1+b_1*height[j]))/(np.sum(height[j:k]) - (k-j)*height[j])
                         
                         a_2 = np.sum(np.multiply(height[j:k], theta[j:k]))/np.sum(height[j:k]) - b_2*np.sum(height[j:k]**2)/np.sum(height[j:k])

                         print np.sum(height[j:k]), np.sum(np.multiply(height[j:k], theta[j:k])), b_2, np.sum(height[j:k]**2)
                         print np.sum(height[j:k]), (k-j)*height[j]

                         b_3 = (np.sum(theta[k:312]) - (312-k)*(a_2+b_2*height[k]))/(np.sum(height[k:312]) - (312-k)*height[k])
                         a_3 = np.sum(np.multiply(height[k:312], theta[k:312]))/np.sum(height[k:312]) - b_2*np.sum(height[k:312]**2)/np.sum(height[k:312])
                         
                         RSS[j, k] = np.sum(np.add(theta[2:j], -(a_1+ b_1*height[2:j]))**2) + np.sum(np.add(theta[j:k], -(a_2+ b_2*height[j:k]))**2) + np.sum(np.add(theta[k:312], -(a_2+ b_2*height[k:312]))**2) 
                         print j, k, a_1, b_1, a_2, b_2, a_3, b_3, RSS[j, k]
     
     [J, K] = np.unravel_index(RSS.argmin(), RSS.shape)
     
     fitvals[:J] = b_1*height[:J] + a_1
     fitvals[J:K] = b_2*height[J:K] + a_2
     fitvals[K:312] = b_3*height[K:312] + a_3


     return fitvals, RSS, J, K                                                                  
                                                                                    
                                                                              
                                                                                       
                                                                                       


#set up plot
theFig = plt.figure(1)
theFig.clf()
theAx = theFig.add_subplot(121)
theAx.set_title('')
theAx.set_xlabel('dtheta/dz (K/m)')
theAx.set_ylabel('z (m)')

theAx1 = theFig.add_subplot(122)
theAx1.set_title('')
theAx1.set_xlabel('')
theAx1.set_ylabel('')

theFig = plt.figure(2)
theFig.clf()
theAx2 = theFig.add_subplot(111)
theAx2.set_title('')
theAx2.set_xlabel('')
theAx2.set_ylabel('')

dump_time_list, time_hrs = Make_Timelists(1, 600, 28800)

dump_time = dump_time_list[29]

for k in range(1):
     
     [wvels, theta, tracer, height] = nc.Get_Var_Arrays("/tera2/nchaparr/Dec252013/runs/sam_case", "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_", dump_time, k+1)

     
     [dvardz, grad_peaks] = nc.Domain_Grad(theta, height) 
     tops_indices=np.where(np.abs(grad_peaks - 0)>0)
     
     #for i in range(tops_indices[0].shape[0]):
     for i in range(1):
          top_index = [tops_indices[0][i], tops_indices[1][i]]
          [i, j] = top_index
          thetavals = theta[:, i, j]
          fitvals, RSS, J, K = get_RSS(thetavals, height) 
          theAx1.plot(thetavals, height[:], 'wo')
          theAx.plot(fitvals, height[:], 'r-')
          
x = np.arange(0, 312, 1)
y = np.arange(0, 312, 1)
X,Y = np.meshgrid(x, y)

     #plt.legend(loc = 'upper center', prop={'size':8})
     #theAx.text(304, 1700, 'h',  fontdict=None, withdash=True)
im = theAx2.pcolor(X, Y, RSS, cmap=cm.bone)
     #im = theAx.pcolor(X, Y, wslice, cmap=cm.bone, vmax=4, vmin=-3) #
bar = plt.colorbar(im)
     #bar.locator = ticker.FixedLocator(tick_list)
     #bar.formatter= ticker.FixedFormatter(label_list)
     #bar.update_ticks()
     #plt.ylim(0, 3200)
     #plt.xlim(0, 4800)

theAx1.set_xlim(300, 310)
theAx1.set_ylim(0, 4000)
theAx.set_ylim(0, 4000)
plt.show()




    
    
