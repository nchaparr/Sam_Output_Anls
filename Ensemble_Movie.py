from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib import cm
import matplotlib.pyplot as plt
#import site
#site.addsitedir('/tera/phil/nchaparr/SAM2/sam_main/python')
#from Percentiles import *
from matplotlib.patches import Patch
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
#import nchap_fun as nc
import matplotlib.animation as animation
#from Ens_Profs import Main_Fun
from Make_Timelist import *
from nchap_class import *

"""
    Profiles/2d ims  at a point
    for a movie
    may be pointless now -- data at longer delta ts
"""


#set up plot

dump_time_list, Times_hrs = Make_Timelists(1, 60, 28800)
dump_time = dump_time_list[59]
i=1


ims = []
ims1 = []
theFig = plt.figure(1)
#theFig1.clf()
#theAx = theFig.add_subplot(111)
#theAx1 = theFig.add_subplot(111)
#theAx.set_title('')
#theAx.set_xlabel('')
#theAx.set_ylabel('')
plt.ylim(0, 1200)
plt.xlim(0, 4790)
#for dump_time in dump_time_list:
for i in range(100):
    dump_time=dump_time_list[i]
    date = "Aug122014" #TODO: this should be an argument passed to Main_Fun
     #pulls data using class Get_Var_Arrays1     
    Vars_File = "/newtera/tera/phil/nchaparr/tera2_cp/nchaparr/"+date+"/runs/sam_case1/OUT_3D/NCHAPP1_testing_doscamiopdata_24_"+ dump_time+".nc"
    #Vars = Get_Var_Arrays1("/newtera/tera/phil/nchaparr/tera2_cp/nchaparr/"+date+"/runs/sam_case", "/OUT_3D/NCHAPP1_testing_doscamiopdata_24_", dump_time)
    #thetapert_list = Vars.get_thetaperts()
    #thetaperts = thetapert_list[0]
    Vars = Dataset(Vars_File, 'r')
    wvel = np.squeeze(Vars.variables['W'][...])
    temp = np.squeeze(Vars.variables['TABS'][...])
    press = np.squeeze(Vars.variables['p'][...])
    height = np.squeeze(Vars.variables['z'][...])
    [znum, ynum, xnum] = temp.shape
    theta = np.zeros_like(temp) #TODO: make this into a function
    thetafact = np.array([(1.0*1000/k)**(1.0*287/1004) for k in press])
    for j in range(znum):
        theta[j, :, :] = temp[j, :, :]*thetafact[j]
    Vars.close()
    print temp.shape
    
    #getting horizontally averaged, ensemble averaged tracer
    #[tracer, theta, height] = Main_Fun(dump_time)        

    #[grad_tracer, tracer_peaks] = nc.Domain_Grad(tracer, height)
    [yindex, xindex] = [13, 44] 
    #print yindex, xindex, tracer_peaks[yindex, xindex]
    i=i+1

    x = np.arange(0, 4800, 25)
    y = height[0:312]
    X,Y = np.meshgrid(x, y)
    dthetadh = np.zeros([znum-1, ynum, xnum])
    for i in range(znum-1):
         dtheta = theta[i+1, :, :] - theta[i, :, :]
         dh = height[i+1] - height[i]
         #print dh, i
         dthetadh[i, :, :] = 1.0*dtheta/dh
          
    #tslice = tracer[0:64, 13, :]
    thetaslice = dthetadh[0:312, 13, :]
    thetaslice = np.ma.masked_where(thetaslice==0, thetaslice)
    thetaslice = np.log(thetaslice)
    
    v_max, v_min = np.amax(thetaslice), np.amin(thetaslice)
    abs_min = np.amin(np.abs(thetaslice))
    print abs_min
    
    #thetaslice = dthetadh[0:312, 13, :]
    #thetaslice = np.ma.masked_inside(thetaslice, -.21, .21)
    im = plt.pcolor(X, Y, thetaslice, cmap=cm.Greys) #, vmax=.01, vmin=.002
    #im = plt.imshow(thetaslice)
    ims.append([im])
    #ims.append(plt.pcolor(X, Y, tslice, norm=plt.Normalize(0, 30)),))
    #ims.append((plt.pcolor(X, Y, tslice, norm=plt.Normalize(0, 30)),))
    #ims.append(plt.pcolor(X, Y, thetaslice, norm=plt.Normalize(0, 30)))
    #ims.append(plt.plot(tracer[:, yindex, xindex], height, 'ko'))
    #ims.append(plt.plot(theta[:210, yindex, xindex], height[:210], 'k-'))    
    #plt.savefig('/tera/phil/nchaparr/python/Plotting/July92013/pngs/for_point_movie/Point_Tracer_'+ str(i)+'.png', bbox_inches=0)

im_ani = animation.ArtistAnimation(theFig, ims, interval=100, repeat_delay=300, blit=True)

#im_ani = animation.ArtistAnimation(theFig, ims, interval=1000, repeat_delay=30000, blit=True)
#im_ani.save('/tera/phil/nchaparr/python/Plotting/July92013/pngs/for_point_movie/im.mp4')

plt.show()

    
    
