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


#set up plot
theFig = plt.figure(1)
theFig.clf()
theAx = theFig.add_subplot(121)
theAx.set_title('')
theAx.set_xlabel('dtheta/dz (K/m)')
theAx.set_ylabel('z (m)')

theAx1 = theFig.add_subplot(122)
theAx1.set_title('')
theAx1.set_xlabel('Theta (K)')
theAx1.set_ylabel('')

dump_time_list, time_hrs = Make_Timelists(1, 600, 28800)


dump_time = dump_time_list[29]

for k in range(10):
     
     [wvels, theta, tracer, height] = nc.Get_Var_Arrays("/tera2/nchaparr/Dec252013/runs/sam_case", "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_", dump_time, k+1)

     #[fit_array, dvardz, grad_peaks] = get_fit(theta, height)
     [dvardz, grad_peaks] = nc.Domain_Grad(theta, height) 
     v_max, v_min, mean, stddev = np.amax(grad_peaks), np.amin(grad_peaks), np.mean(grad_peaks), np.std(grad_peaks)
     #label_list, tick_list = get_ticks(mean, stddev,v_max, v_min)
     print 'max min std', v_max, v_min, mean, stddev
     #get slice of w, at domain averaged h
     yavh = np.mean(grad_peaks, axis = 0)
     avh = np.mean(yavh, axis = 0)
     slice_index = np.where(height > 1500)[0]

     #tops = tracer_peaks[np.abs(tracer_peaks - 1570)<10]
     #print tops.shape, tops
     tops_indices=np.where(np.abs(grad_peaks - 0)>0)
     
     for i in range(tops_indices[0].shape[0]):
     #for i in range(1):
          top_index = [tops_indices[0][i], tops_indices[1][i]]
          [i, j] = top_index
          if np.mod(i, 10) ==0:
               print i, j
          #print dvardz[50-1:0:-1, i, j].shape, dvardz[:,i, j].shape, dvardz[-1:50:-1, i, j].shape
          #s = np.r_[dvardz[50-1:0:-1, i, j], dvardz[:,i, j], dvardz[-1:50:-1, i, j]]
          #print s.shape
          #w = np.ones(50, 'd')
          #y=np.convolve(w/w.sum(), s, mode='valid')
          #print y.shape
          
          #gradmax = np.amax(dvardz[: , i, j])
          #index = np.where(dvardz[:, i, j] - gradmax == 0)[0][0]
          #index2 = np.argmax(dvardz[:, i, j])
          #print index, index2, height[index], gradmax, grad_peaks[i,j]
     #top_index = [24, 68]

          #print top_index, top_index[0]*25, top_index[1]*25, tracer_peaks[top_index[0], top_index[1]], v_max
     #theAx.plot(np.zeros(312)+302, height, 'k--')

     #theAx.plot(wvels[:, top_index[0], top_index[1]]+302, height, 'r--', label='w')
          thetavals = theta[:, top_index[0], top_index[1]]
          #s1 = np.r_[theta[10-1:0:-1, i, j], theta[:,i, j], theta[-1:10:-1, i, j]]
          #print s1.shape
          #w1 = np.ones(10, 'd')
          #y1=np.convolve(w1/w1.sum(), s1, mode='valid')
          #print y.shape
          
          #fitvals = fit_array[:, top_index[0], top_index[1]]
          gradvals = dvardz[:, top_index[0], top_index[1]]
          #dheight = np.gradient(height)
          #dvar = np.gradient(thetavals)
          #dvardz = np.divide(dvar, dheight)
          #print ddvardzdz.shape
          #stepvals = steps[:, top_index[0], top_index[1]]
          theAx1.plot(thetavals, height[:], 'wo')
          #theAx1.plot(y1[:312], height[:], 'b-')
          gradref0 = np.zeros_like(gradvals)
          gradref1 = np.zeros_like(gradvals) + .0025
          theAx.plot(gradvals, height[:], 'r-')
          #theAx.plot(wvels[:, top_index[0], top_index[1]], height[:], 'b-')
          theAx.plot(gradref0, height[:], 'k--')
          theAx.plot(gradref1, height[:], 'k--')          
          #theAx.plot(dvardz, height[:], 'k-')
          #theAx.plot(y[:312], height[:], 'b-')
          #theAx1.plot(fitvals+.5, height[:], 'b-')
          
     #filler_array = np.zeros([64, 192])
     #wslice = wvels[slice_index[0], :, :]
     #wslice = np.vstack((wslice, filler_array))
     #n, bins, patches = theAx.hist(tracer_peaks, bins=20)
     #height_bin_vols = nc.Bin_Peaks(grad_peaks, height)
     #theAx.bar(height, height_bin_vols)
     #grad_peaks = np.vstack((grad_peaks, filler_array))
     #x = np.arange(0, 4800, 25)
     #y = np.arange(0, 4800, 25)
     #X,Y = np.meshgrid(x, y)

     #plt.legend(loc = 'upper center', prop={'size':8})
     #theAx.text(304, 1700, 'h',  fontdict=None, withdash=True)
     #im = theAx.pcolor(X, Y, grad_peaks, cmap=cm.bone, vmax=v_max, vmin=v_min)
     #im = theAx.pcolor(X, Y, wslice, cmap=cm.bone, vmax=4, vmin=-3) #
     #bar = plt.colorbar(im)
     #bar.locator = ticker.FixedLocator(tick_list)
     #bar.formatter= ticker.FixedFormatter(label_list)
     #bar.update_ticks()
     #plt.ylim(0, 3200)
     #plt.xlim(0, 4800)

theAx1.set_xlim(300, 310)
theAx1.set_ylim(0, 4000)
theAx.set_ylim(0, 4000)
plt.show()




    
    
