import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from matplotlib import rcParams
rcParams.update({'font.size': 10})

"""

   For plotting the (scaled) temperature gradient and flux profiles.

"""

date = "Mar12014"
sfc_flx = 60
gamma = .01

Fig1 = plt.figure(1)
Fig1.clf()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


Ax = Fig1.add_subplot(131)
#Ax.set_title( r'$\theta$', fontsize=20)
#Ax.set_title( r'$\frac{\partial \theta}{\partial z}$', fontsize=20)
#Ax.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
Ax.set_xlabel(r"$\overline{\theta}$", fontsize=20)
Ax.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#Ax.set_ylabel(r"$z$", fontsize=20)
plt.xlim(300, 312)
#plt.ylim(100, 1600)
plt.ylim(0.1, 1.4)


Ax1 = Fig1.add_subplot(132)
#Ax1.set_title( r'$Scaled \ \frac{\partial \theta}{\partial z}$', fontsize=20)
#Ax1.set_title( r'$\frac{\partial \theta}{\partial z}$', fontsize=20)
Ax1.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
#Ax1.set_xlabel(r"$\frac{\partial \theta}{\partial z}$", fontsize=20)
#Ax1.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#start, end = -.025, .025
start, end = -1, 2.5
Ax1.set_xticks([.02, 1])
#Ax1.set_ylabel(r"$z$", fontsize=20)
#plt.xlim(-.025, .025)
#plt.xlim(-1, 2.5)
#plt.ylim(100, 1600)
plt.ylim(0.1, 1.4)

Ax2 = Fig1.add_subplot(133)
#Ax2.set_title(r"$\overline{w^{'} \theta^{'}}$", fontsize=20)
#Ax2.set_title(r"$Scaled \ \overline{w^{'} \theta^{'}}$", fontsize=20)
#Ax2.set_xlabel(r"$\overline{w^{'}\theta^{'}}$", fontsize=20)
Ax2.set_xlabel(r"$\frac{\overline{w^{'}\theta^{'}}}{\overline{w^{'}\theta^{'}}_{0}}$", fontsize=20)
#start, end = -.06, .14
#start, end = -.4, 1.2
#Ax2.set_xticks(np.arange(start, end, 1.0*(end-start)/5))

#Ax2.set_ylabel(r"$z$", fontsize=20)
#Ax2.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#plt.ylim(100, 1600)
#plt.xlim(-.06, .14)
#plt.xlim(-.4, 1.2)
plt.ylim(0.1, 1.4)
dump_time_list, Times = Make_Timelists(1, 600, 28800)
 
theta_file_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/theta_bar"+ dump_time for dump_time in dump_time_list]
press_file_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/press"+ dump_time for dump_time in dump_time_list]
flux_file_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/wvelthetapert"+ dump_time for dump_time in dump_time_list]
height_file = "/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/heights0000000600"
AvProfVars = np.genfromtxt("/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/AvProfLims")

#loop over text files files
for i in range(len(theta_file_list)):
    
    theta = np.genfromtxt(theta_file_list[i])
    height = np.genfromtxt(height_file)
    press = np.genfromtxt(press_file_list[i])
    rhow = nc.calc_rhow(press, height, theta[0])
    wvelthetapert = np.genfromtxt(flux_file_list[i])
    wvelthetapert[0] = np.nan
    #Now for the gradients
    dheight = np.diff(height)
    dtheta = np.diff(theta)      
    dthetadz = np.divide(dtheta, dheight)
           
    element0 = np.array([0])
    dthetadz=np.hstack((dthetadz, element0))
        
    #only need up to 2500meters
    top_index = np.where(abs(2545 - height) < 40.)[0][0]

    #where gradient is max, and flux is min
    #print AvProfVars[:,1].shape, height.shape
    scaled_height = [1.0*h/AvProfVars[i,1] for h in height]

    fluxes = np.multiply(wvelthetapert, rhow)*1004.0/sfc_flx
    
    if np.mod(i+1, 1) == 0:
    #if i > 14 and i < 21:
        
        fluxes[0] = np.nan
        zeros = np.zeros_like(height)

        Ax.plot(theta, scaled_height, '-') #, label = str(Times[i])+'hrs'
        
        Ax1.plot(1.0*dthetadz/gamma, scaled_height, '-', label = str(Times[i])+'hrs')
        
        Ax2.plot(fluxes, scaled_height, '-', label = str(Times[i])+'hrs')    
    
array = np.genfromtxt('/newtera/tera/phil/nchaparr/python/Pert_Files/snd')
    
height_0 = array[:, 0]
theta_0 = array[:, 1]
f=interp1d(height_0, theta_0)

#Now plot inital sounding
top_index = np.where(height <= 2500)[0][-1]
theta_0 = f(height[0:top_index])
dtheta0 = np.diff(theta_0)
dthetadz0 = np.divide(dtheta0, dheight[0:top_index-1])
element0 = np.array([.005])
dthetadz0=np.hstack((element0, dthetadz0))

#Ax1.plot(dthetadz0, scaled_height[0:top_index], '--', label = 'Initial Sounding')
Ax1.plot(zeros, height)#zeros line for reference
#Ax1.plot(gamma, scaled_height)#zeros line for reference
Ax1.plot(zeros+1, height, 'k-')#zeros line for reference
Ax1.plot(zeros+.02, height, 'k-')#zeros line for reference
#Ax2.plot(zeros, height)#zeros line for reference
Ax2.plot(zeros, height)#zeros line for reference 
#plt.legend(loc = 'Lower right', prop={'size':8})
#Ax2.plot(theta_0, scaled_xheight[0:top_index], '--', label = 'Initial Sounding')#"
#plt.xlim(300, 310)
#plt.legend(loc = 'upper right', prop={'size':8})
plt.show()
#Fig1.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/theta_flux_profs.png")





    
    
