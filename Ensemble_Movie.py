from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import nchap_fun as nc
import matplotlib.animation as animation
from Make_Timelist import *
from nchap_class import *
import fastfit as fsft
import datetime


"""
    Profiles/2d ims  at a point/xz-slice
    for a movie
    for defense
"""


def Get_CBLHeights_thetas(heights, press, thetas, gamma, top_index):
    """
    Gets heights based on dthetdz
    
    Arguments:
     -- 
     --

    Returns:
     --   

    """
    
    dheight = np.diff(heights)
    dtheta = np.diff(thetas)
    dthetadz=np.divide(dtheta, dheight)
    element0=np.array([0])
    dthetadz=np.hstack((element0, dthetadz))*1.0/gamma
    rhow=calc_rhow(press, heights, thetas[0])

    #where gradient is greater than zero    
    for j in range(len(dthetadz[:top_index])-1):
        if (dthetadz[j+1] >.03) and (dthetadz[j] >= 0):
            dtheta_index_b = j+1
            break

    #where gradient resumes as gamma
    dtheta_index_t = 999
    for k in range(len(dthetadz[:top_index])-1):
        #print dthetadz[k-1], dthetadz[k+1], dthetadz[k+2]
        #print ""
        #print np.abs(dthetadz[k+1]-1), np.abs(dthetadz[k+2]-1)
        if np.abs(dthetadz[k+2]-1)<.03 and np.abs(dthetadz[k+1]-1)<.03 and dthetadz[k-1]>1:            
            dtheta_index_t = k+1                        
            break

    #Hacky fix for when the upper theta gradient profiles are wonky
    if dtheta_index_t == 999:
         for k in range(len(dthetadz[:top_index])-1):
           #   print dthetadz[k-1], dthetadz[k+1], dthetadz[k+2]
           #   print ""
           #   print np.abs(dthetadz[k+1]-1), np.abs(dthetadz[k+2]-1)
              if np.abs(dthetadz[k+2]-1)<.04 and np.abs(dthetadz[k+1]-1)<.04 and dthetadz[k-1]>1:            
                   dtheta_index_t = k+1                        
                   break

    eltop_dthetadz = heights[dtheta_index_t]
    elbot_dthetadz = heights[dtheta_index_b]
    
    h = heights[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]]
    
    deltatheta = thetas[dtheta_index_t] - thetas[dtheta_index_b]
    mltheta = np.mean(thetas[0:dtheta_index_b])
    
    return [elbot_dthetadz, h, eltop_dthetadz, deltatheta, mltheta]

def get_dthetadh(theta, height):
    dthetadh = np.zeros([znum-1, ynum, xnum])
    for i in range(znum-1):
        dtheta = theta[i+1, :, :] - theta[i, :, :]
        dh = height[i+1] - height[i]
            #print dh, i
        dthetadh[i, :, :] = 1.0*dtheta/dh
        return dthetadh

def get_ML_Heights(thetaslice, height):
        ML_Heights = np.empty([120])
        for i in range(120): #change to 128          
             #for j in range(192):#change to 192
             #top = np.where(np.abs(height-2300)<100)[0][0]
            RSS, J, K = fsft.get_fit(thetaslice[:, i], height, 240)
            ML_Heights[i] = height[J]
            #print RSS, J, K, height[J]
            return ML_Heights
    
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
plt.xlim(0, 3000)
#for dump_time in dump_time_list:
for i in range(5):
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
    
    
    #getting horizontally averaged average profile
    #[tracer, theta, height] = Main_Fun(dump_time)
    av_theta = nc.Horizontal_Average(theta)
    [elbot_dthetadz, h, eltop_dthetadz, deltatheta, mltheta] = Get_CBLHeights_thetas(height, press, av_theta, .01, 241)

    #[grad_tracer, tracer_peaks] = nc.Domain_Grad(tracer, height)
    [yindex, xindex] = [13, 44] 
    #print yindex, xindex, tracer_peaks[yindex, xindex]
    i=i+1

    x = np.arange(0, 3000, 25)
    y = height[0:312]
    X,Y = np.meshgrid(x, y)

    dthetadh = get_dthetadh(theta, height) 
    #tslice = tracer[0:64, 13, :]
    #thetaslice = dthetadh[0:312, 13, :]
    #thetaslice = np.ma.masked_where(thetaslice==0, thetaslice)
    
    #thetaslice = np.log(thetaslice)
    
    #v_max, v_min = np.amax(thetaslice), np.amin(thetaslice)
    #abs_min = np.amin(np.abs(thetaslice))
    #print abs_min
    
    thetaslice = theta[0:312, 114, 40:160]
    t0 = datetime.time(1, 2, 3)
    print t0.second
    t0_sec = t0.second
    
    ML_Heights = get_ML_Heights(thetaslice, height)
    
    tnow=datetime.time(1, 2, 3)
    tnow_sec = tnow.second
    print tnow_sec - t0_sec 
    
    thetaslice = np.ma.masked_where(thetaslice>mltheta+2, thetaslice)

    print thetaslice.shape
    #thetaslice = np.ma.masked_inside(thetaslice, mltheta-5, mltheta+5)
    #im = plt.pcolor(X, Y, thetaslice-mltheta, cmap=cm.bone) #, vmax=.01, vmin=.002
    #im = plt.imshow(thetaslice)
    #ims.append([im])
    #ims.append(plt.pcolor(X, Y, tslice, norm=plt.Normalize(0, 30)),))
    #ims.append((plt.pcolor(X, Y, tslice, norm=plt.Normalize(0, 30)),))
    #ims.append(plt.pcolor(X, Y, thetaslice, norm=plt.Normalize(0, 30)))
    ims.append(plt.plot(x, ML_Heights, 'k-'))
    #ims.append(plt.plot(theta[:210, yindex, xindex], height[:210], 'k-'))    
    #plt.savefig('/tera/phil/nchaparr/python/Plotting/July92013/pngs/for_point_movie/Point_Tracer_'+ str(i)+'.png', bbox_inches=0)

im_ani = animation.ArtistAnimation(theFig, ims, interval=300, repeat_delay=1000, blit=True)

#im_ani = animation.ArtistAnimation(theFig, ims, interval=1000, repeat_delay=30000, blit=True)
mywriter = animation.FFMpegWriter()
im_ani.save('im.mp4', writer=mywriter)

plt.show()

    
    
