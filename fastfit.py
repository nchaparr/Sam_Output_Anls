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
     

     return RSS, j, k                                                              
                                                                               


    
    
