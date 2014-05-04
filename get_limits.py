import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from nchap_class import *
from matplotlib import rcParams
rcParams.update({'font.size': 10})



"""calculates temperature gradients (discrete) from txt files inturn from ensemble run 3D files
   gets levels where gradient exceeds zero, and where it resumes gamma, and zero crossings for fluxes
   and the maximal points withing the entrainment region.
   Gets Delta Theta and the Mixed Layer average Theta.
   Dumps them in a text file.
   Calcs and dumps rino, invrino, wstar
"""
#to be changed for each run
rundate = 'Mar52014'
gamma = .01
flux_s = 150

#output times
dump_time_list, Times = Make_Timelists(1, 600, 28800)
Times = np.array(Times)  

files = For_Plots(rundate)

theta_file_list = [files.get_file(dump_time, "theta_bar") for dump_time in dump_time_list]
press_file_list = [files.get_file(dump_time, "press") for dump_time in dump_time_list]
#press_file_list = ["/tera/phil/nchaparr/python/Plotting/" + rundate + "/data/press"+ dump_time for dump_time in dump_time_list]
flux_file_list = [files.get_file(dump_time, "wvelthetapert") for dump_time in dump_time_list]
#flux_file_list = ["/tera/phil/nchaparr/python/Plotting/" + rundate + "/data/wvelthetapert"+ dump_time for dump_time in dump_time_list]
#height_file = "/tera/phil/nchaparr/python/Plotting/" + rundate + "/data/heights0000000600"
height_file = files.get_file("0000000600", "heights")

AvProfLims = []
invrinos = []
#loop over text files files
for i in range(len(theta_file_list)):
    print i, theta_file_list[i]
    theta = np.genfromtxt(theta_file_list[i])
    print theta.shape
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
    top_index = np.where(abs(1700 - height) < 26.)[0][0] #may need to be higher (e.g. for 60/2.5)
       
    #TODO: add method to For_Plots, get dtheta heights, get flux heights
    #where gradient is greater than zero    
    #for j in range(len(dthetadz[:top_index])-1):
    #    if (dthetadz[j+1] >.0002) and (dthetadz[j] >= 0):
    #        dtheta_index_b = j+1
    #        break

    #where gradient resumes as gamma    
    #for k in range(len(dthetadz[:top_index])-1):        
    #    if np.abs(dthetadz[k+2]-gamma)<.0002 and np.abs(dthetadz[k+1]-gamma)<.0002 and dthetadz[k-1]>gamma:            
    #        dtheta_index_t = k+1                        
    #        break
    
    #now fluxes    
    #fluxes = np.multiply(wvelthetapert, rhow)*1004.0
    #for l in range(len(dthetadz)-1):
    #    if (fluxes[l+1] <= .0) and (fluxes[l] > 0):
    #        flux_index_b = l+1
    #        break
        
    #for m in range(len(dthetadz[0:top_index])-1):
    #    if (abs(fluxes[m+1]) < 0.5) and (fluxes[m] < 0) and (fluxes[m-1] < 0):
    #        flux_index_t = m+1
    #        break
        
    #eltop_dthetadz = height[dtheta_index_t]
    #elbot_dthetadz = height[dtheta_index_b]
    #eltop_flux = height[flux_index_t]
    #elbot_flux = height[flux_index_b]
    #h = height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]]
    #h_flux = height[np.where(wvelthetapert - np.amin(wvelthetapert) == 0)[0][0]]

    print height.shape, press.shape, theta.shape, wvelthetapert.shape, gamma, top_index
    
    [elbot_dthetadz, h, eltop_dthetadz, elbot_flux ,h_flux  ,eltop_flux, deltatheta, mltheta]= nc.Get_CBLHeights(height, press, theta, wvelthetapert, gamma, top_index)
    print elbot_dthetadz, h, eltop_dthetadz, elbot_flux ,h_flux  ,eltop_flux, deltatheta, mltheta
    #mltheta = np.mean(theta[0:elbot_dthetadz])
    #print theta[eltop_dthetadz], elbot_dthetadz
    #deltatheta = theta[eltop_dthetadz] - theta[elbot_dthetadz]
    
    delta_h=eltop_dthetadz - elbot_dthetadz
    #TODO: this can be tidied up, ie name valriables and pass the named variables to calc_rino    
    #print i, height[dtheta_index_b], height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]], height[dtheta_index_t], height[flux_index_b], height[np.where(wvelthetapert - np.amin(wvelthetapert) == 0)[0][0]], height[flux_index_t]
    #print i, height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]], np.mean(theta[0:dtheta_index_b]), 1.0*flux_s/(rhow[0]*1004), -theta[dtheta_index_b]+theta[dtheta_index_t]
    
    [rino, invrino, wstar, S, pi3, pi4] =  nc.calc_rino(h, mltheta, 1.0*flux_s/(rhow[0]*1004), deltatheta, gamma, delta_h)

    AvProfLims.append([elbot_dthetadz, h, eltop_dthetadz, elbot_flux, h_flux, eltop_flux, deltatheta, mltheta])
    tau = 1.0*h/wstar
    invrinos.append([rino, invrino, wstar, S, tau, mltheta, deltatheta, pi3, pi4])

files.save_file(np.array(AvProfLims), "AvProfLims")
files.save_file(np.array(invrinos), "invrinos")

#np.savetxt('/tera/phil/nchaparr/python/Plotting/' + rundate + '/data/AvProfLims', np.array(AvProfLims), delimiter=' ')
#np.savetxt('/tera/phil/nchaparr/python/Plotting/' + rundate + '/data/invrinos', np.array(invrinos), delimiter=' ')




    
    
