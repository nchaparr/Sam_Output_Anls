import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
import scipy.special as scsp
from Make_Timelist import *

#import site
#site.addsitedir('/tera/phil/nchaparr/python')
import nchap_fun as nc


#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
#import nchap_fun as nc

"""
   copied from get_rino
   for looking at the root mean w squared behaviour as a measure of TKE
   
"""
F0 = 75 #j/sm2
rho = 1 #Kg/m3
cp = 1004 #j/KgK
date = "Mar52014"
gamma = .01

#TODO: See if plots can be set up using a function

Fig1 = plt.figure(1)
Fig1.clf()

Ax1 = Fig1.add_subplot(121)
Ax1.set_title('')
Ax1.set_ylabel(r'$\frac{z}{h}$', fontsize=20)
Ax1.set_xlabel(r'$\sqrt[]{^{u\prime 2}} \ (m/s)$', fontsize=20)
plt.ylim(.1, 1.4)

Ax2 = Fig1.add_subplot(122)
Ax2.set_xlabel(r'$\frac{\sqrt[]{u^{\prime 2}}}{w^{*}}$', fontsize=20)

#plt.xlim(0, 5)
plt.ylim(.1, 1.4)

#plt.xlim(0, 5)


#create lists of txt file to loop over

dump_time_list, Times = Make_Timelists(1, 600, 28800)

theta_file_list = ["/tera/phil/nchaparr/python/Plotting/"+date+"/data/theta_bar"+ dump_time for dump_time in dump_time_list]
flux_file_list = ["/tera/phil/nchaparr/python/Plotting/"+date+"/data/wvelthetapert"+ dump_time for dump_time in dump_time_list]
height_file = "/tera/phil/nchaparr/python/Plotting/"+date+"/data/heights0000000600"
press_file_list = ["/tera/phil/nchaparr/python/Plotting/"+date+"/data/press"+ dump_time for dump_time in dump_time_list]
rmwsq_file_list = ["/tera/phil/nchaparr/python/Plotting/"+date+"/data/rootmeanusq"+ dump_time for dump_time in dump_time_list]



BLheights = []
BLheights0=[]
BLheights1=[]
BLheights10=[]
BLheights11=[]
invrinos = []
#BLheightscheck=[]
theta_jump=[]
theta_ml=[]
wstars = []
maxrmwsqs = []
rmwsqs_at_h = []
wthetas_s = []

#loop over text files files
for i in range(len(theta_file_list)):
    theta = np.genfromtxt(theta_file_list[i])
    height = np.genfromtxt(height_file)
    press = np.genfromtxt(press_file_list[i])
    rmwsq = np.genfromtxt(rmwsq_file_list[i])
    #print press.shape
    rhow = nc.calc_rhow(press, height, theta[0])
    #Now for the gradients
    dheight = np.diff(height)
    dtheta = np.diff(theta)
    
    dthetadz = np.divide(dtheta, dheight)
    
    element0 = np.array([0])
    dthetadz=np.hstack((element0, dthetadz))
        
    #only need up to 2500meters
    top_index = np.where(abs(1500 - height) < 50.)[0][0]

    #TODO: see test_lambda
    max_dtheta = np.max(dthetadz[0:top_index]) 
    max_dtheta_index = np.where(abs(dthetadz - max_dtheta)<.00006)[0][0]

    maxrmwsq = np.max(rmwsq[0:top_index])
    rmwsq_at_h = rmwsq[max_dtheta_index]
    #print max_dtheta_index
    BLheights.append(height[max_dtheta_index])
    #where gradient is greater than zero    
    for j in range(len(dthetadz)-1):
         if (dthetadz[j+1] >.0002) and (dthetadz[j] >= 0):
            dtheta_index_b = j+1
            break

    #where gradient resumes as gamma    
    for k in range(len(dthetadz)-1):
         if np.abs(dthetadz[k+1]-gamma)<.0002 and np.abs(dthetadz[k]-gamma)<.0002 and np.abs(dthetadz[k]-gamma)<.0002:
            dtheta_index_t = k+1
            break

    #now fluxes
    wvelthetapert = np.genfromtxt(flux_file_list[i])
    #fluxes = wvelthetapert
    fluxes = np.multiply(wvelthetapert, rhow)*1004 #commented this out for comp to test_lambda    
    for l in range(len(dthetadz)-1):
        if (fluxes[l+1] <= .0) and (fluxes[l] > 0):
            flux_index_b = l+1
            break

    #print np.where(abs(.00 - fluxes) < .3)[0][0:3]
    for m in range(len(dthetadz)-1):
        if (abs(fluxes[m+1]) < 0.8) and (fluxes[m] < 0):
            flux_index_t = m+1
            break
    
    BLheights0.append(height[dtheta_index_b])
    BLheights1.append(height[dtheta_index_t])    

    BLheights10.append(height[flux_index_b])
    BLheights11.append(height[flux_index_t])
           
    theta_jump.append(theta[dtheta_index_t]-theta[dtheta_index_b])
    theta_ml.append(np.mean(theta[0:dtheta_index_b]))    
       
    #plt.legend(loc = 'upper right', prop={'size':8})
    #wtheta_s = 1.0*60*scsp.erf(Times[i]/(2.5*np.sqrt(2)))/(1004*rhow[0])
    wtheta_s = 150/(1004*rhow[0])
    #print 'rhow', rhow[0], 1.0*60*scsp.erf(Times[i]/(2.5*np.sqrt(2)))
    print BLheights[i], theta_ml[i], wtheta_s, theta_jump[i]
    [rino, invrino, wstar, S, pi3, pi4] =  nc.calc_rino(BLheights[i], theta_ml[i], wtheta_s, theta_jump[i], gamma, height[dtheta_index_t]-height[dtheta_index_b])
    invrinos.append(invrino)
    wstars.append(wstar)
    wthetas_s.append(wtheta_s)
    maxrmwsqs.append(maxrmwsq)
    rmwsqs_at_h.append(rmwsq_at_h)

    if np.mod(i+1, 6)==0:
        #print maxrmwsq, wstar, 1.0*maxrmwsq/wstar
        Ax1.plot(rmwsq, 1.0*height/BLheights[i],label = str(Times[i]) + 'hrs')
        Ax2.plot(1.0*rmwsq/wstar, 1.0*height/BLheights[i], label = str(Times[i]) + 'hrs')
    #plt.legend(loc = 'upper right', prop={'size':8})

        
    
    #TODO: Need to make sure sfc flux is ok (w'theta'0)

#zeros = np.zeros_like(height[0:top_index])

#dhdt = nc.get_dhdt(np.array(BLheights), np.array(Times))
#we = np.divide(dhdt, wstars[0:48])
#Ax2.plot(zeros, height[0:top_index])
plt.legend(loc = 'upper right', prop={'size':8})
plt.show()
#Times = [Time*5 for Time in Times]
#plot the heights vs time
#Fig2 = plt.figure(2)
#Fig2.clf()
#Ax3 = Fig2.add_subplot(111)
#deltaH = [1.0*(i - j)/i for i,j  in zip(BLheights11, BLheights10)]
#print deltaH
#Ax3.plot(Times, theta_jump, 'yo', label = 'theta0')
#print len(Times), len(invrinos)
#Ax3.plot(Times, wstars, 'ko')
#Ax3.plot(Times, wthetas_s,'y*')
#Ax3.plot(Times, theta_jump,'k*', label = 'flux1')
#plt.ylim(-600, 600)
#plt.xlim(0, .04)
#plt.legend(loc = 'lower right', prop={'size':8})
#Ax3.set_title('Wstar vs Time')
#Ax3.set_ylabel('Wstar (m/s)')
#Ax3.set_xlabel('Time (hrs)')
#plt.show()



    
    
