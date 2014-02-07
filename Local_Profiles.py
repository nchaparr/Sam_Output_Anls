from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib import cm
import matplotlib.pyplot as plt

import site
site.addsitedir('/tera/phil/nchaparr/python')
#from Percentiles import *
from matplotlib.patches import Patch

#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc

import matplotlib.animation as animation
from Make_Timelist import *

"""
    Profiles at a point
    for a movie
    only usefull with small delta t
"""


#set up plot
dump_time_list, time_hrs = Make_Timelists(1, 60, 28800)
#dump_time = dump_time_list[59]
i=1


ims = []
ims1 = []
theFig = plt.figure()
#theFig1.clf()
#theAx = theFig.add_subplot(111)
#theAx1 = theFig.add_subplot(111)
#theAx.set_title('')
#theAx.set_xlabel('')
#theAx.set_ylabel('')

for dump_time in dump_time_list:
    
    [wvels, theta, tracer, height] = nc.Get_Var_Arrays('/tera/phil/nchaparr/python/Plotting/Sep302013/case', '/OUT_3D/NCHAPP1_testing_doscamiopdata_16_', dump_time, 1)
    [grad_theta, theta_peaks] = nc.Domain_Grad(theta, height)
    [yindex, xindex] = [13, 44] 
    #print yindex, xindex, tracer_peaks[yindex, xindex]
    i=i+1

    x = np.arange(0, 3200, 50)
    y = height[0:70]
    X,Y = np.meshgrid(x, y)
    
    tslice = theta[0:70, 13, :]-273
    
    #ims.append((plt.pcolor(X, Y, tslice),))
    ims.append(plt.plot(theta[:, yindex, xindex], height, 'ko'))
    #plt.savefig('/tera/phil/nchaparr/python/Plotting/July92013/pngs/for_point_movie/Point_Tracer_'+ str(i)+'.png', bbox_inches=0)

im_ani = animation.ArtistAnimation(theFig, ims, interval=500, repeat_delay=3000, blit=True)

#im_ani = animation.ArtistAnimation(theFig, ims, interval=1000, repeat_delay=30000, blit=True)
#im_ani.save('/tera/phil/nchaparr/python/Plotting/July92013/pngs/for_point_movie/im.mp4')

plt.show()

    
    
