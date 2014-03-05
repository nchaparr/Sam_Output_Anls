import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from matplotlib import rcParams
rcParams.update({'font.size': 10})

#TODO: flux and gamm need to be passed as arguments

"""calculates temperature gradients (discrete) from txt files inturn from ensemble run 3D files
   gets levels where gradient exceeds zero, and where it resumes gamma, and zero crossings for fluxes
   and the maximal points withing the entrainment region.
   Gets Delta Theta and the Mixed Layer average Theta.
   Dumps them in a text file.
   Calcs and dumps rino, invrino, wstar
"""
rundate = 'Jan152014_1'
gamma = .005
flux_s = 150
dump_time_list, Times = Make_Timelists(1,600, 28800)
Times = np.array(Times)  
 
theta_file_list = ["/tera/phil/nchaparr/python/Plotting/" + rundate + "/data/theta_bar"+ dump_time for dump_time in dump_time_list]
press_file_list = ["/tera/phil/nchaparr/python/Plotting/" + rundate + "/data/press"+ dump_time for dump_time in dump_time_list]
flux_file_list = ["/tera/phil/nchaparr/python/Plotting/" + rundate + "/data/wvelthetapert"+ dump_time for dump_time in dump_time_list]
height_file = "/tera/phil/nchaparr/python/Plotting/" + rundate + "/data/heights0000000600"

AvProfLims = []
invrinos = []
#loop over text files files
for i in range(len(theta_file_list)):
    
    theta = np.genfromtxt(theta_file_list[i])
    height = np.genfromtxt(height_file)
    press = np.genfromtxt(press_file_list[i])
    rhow = nc.calc_rhow(press, height, theta[0])
    wvelthetapert = np.genfromtxt(flux_file_list[i])

    #Now for the gradients
    dheight = np.diff(height)
    dtheta = np.diff(theta)
    
    dthetadz = np.divide(dtheta, dheight)
    
    element0 = np.array([0])
    dthetadz=np.hstack((element0, dthetadz))
        
    #only need up to 1900meters
    top_index = np.where(abs(1625 - height) < 26.)[0][0]
       
    #TODO: see test_lamda
    #where gradient is greater than zero    
    for j in range(len(dthetadz)-1):
        if (dthetadz[j+1] >.0002) and (dthetadz[j] >= 0):
            dtheta_index_b = j+1
            break

    #where gradient resumes as gamma    
    for k in range(len(dthetadz[:top_index])-1):        
        if np.abs(dthetadz[k+2]-gamma)<.0002 and np.abs(dthetadz[k+1]-gamma)<.0002 and dthetadz[k-1]>gamma:            
            dtheta_index_t = k+1                        
            break
    
    #now fluxes    
    fluxes = np.multiply(wvelthetapert, rhow)*1004.0
    for l in range(len(dthetadz)-1):
        if (fluxes[l+1] <= .0) and (fluxes[l] > 0):
            flux_index_b = l+1
            break
        
    for m in range(len(dthetadz[0:top_index])-1):
        if (abs(fluxes[m+1]) < 0.5) and (fluxes[m] < 0) and (fluxes[m-1] < 0):
            flux_index_t = m+1
            break
        
    print height[dtheta_index_b], height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]], height[dtheta_index_t], height[flux_index_b], height[np.where(wvelthetapert - np.amin(wvelthetapert) == 0)[0][0]], height[flux_index_t]
    print height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]], np.mean(theta[0:dtheta_index_b]), 1.0*flux_s/(rhow[0]*1004), -theta[dtheta_index_b]+theta[dtheta_index_t]
    mltheta = np.mean(theta[0:dtheta_index_b])
    deltatheta = -theta[dtheta_index_b]+theta[dtheta_index_t]
    
    [rino, invrino, wstar, S] =  nc.calc_rino(height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]], np.mean(theta[0:dtheta_index_b]), 1.0*flux_s/(rhow[0]*1004), -theta[dtheta_index_b]+theta[dtheta_index_t], .01)

    AvProfLims.append([height[dtheta_index_b], height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]], height[dtheta_index_t], height[flux_index_b], height[np.where(wvelthetapert - np.amin(wvelthetapert) == 0)[0][0]], height[flux_index_t], -theta[dtheta_index_b]+theta[dtheta_index_t], np.mean(theta[0:dtheta_index_b])])
    tau = 1.0*height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]]/wstar
    invrinos.append([rino, invrino, wstar, S, tau, mltheta, deltatheta])
    
np.savetxt('/tera/phil/nchaparr/python/Plotting/' + rundate + '/data/AvProfLims', np.array(AvProfLims), delimiter=' ')
np.savetxt('/tera/phil/nchaparr/python/Plotting/' + rundate + '/data/invrinos', np.array(invrinos), delimiter=' ')




    
    
