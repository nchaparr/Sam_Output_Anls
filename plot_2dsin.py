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
    To get a visual of a 2d sinewave, to understand 2d fts          
"""

#set up plot
theFig = plt.figure(1)
theFig.clf()
theAx = theFig.add_subplot(111)
theAx.set_title('Sine*Sine  - Cos*Cos: k = 1')
theAx.set_xlabel('m')
theAx.set_ylabel('m')

x = np.arange(0, 3200, 25)
y = np.arange(0, 3200, 25)
X,Y = np.meshgrid(x, y)

print ''
print np.sin(2*np.pi*x/3200).shape, np.sin(2*np.pi*y/3200).shape
print ''
print np.sin(2*np.pi*x/3200)*np.sin(2*np.pi*y/3200).shape
print ''
Sin2D = np.sin(2*np.pi*X/3200)*np.sin(2*np.pi*Y/3200) - np.cos(2*np.pi*X/3200)*np.cos(2*np.pi*Y/3200) 
print Sin2D.shape
im = theAx.pcolor(X, Y, Sin2D, cmap=cm.bone)
bar = plt.colorbar(im)
     #bar.locator = ticker.FixedLocator(tick_list)
     #bar.formatter= ticker.FixedFormatter(label_list)
     #bar.update_ticks()
     #plt.ylim(0, 3200)
     #plt.xlim(0, 4800)

#theAx1.set_xlim(300, 310)
#theAx1.set_ylim(0, 4000)
theAx.set_ylim(0, 3200)
theAx.set_xlim(0, 3200)
plt.show()




    
    
