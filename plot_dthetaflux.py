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

label_list = ['100/10', '100/5', '60/5', '60/2.5', '150/5', '60/10', '150/10']
legend_list = ['kv', 'ko', 'yo', 'y*', 'ro', 'yv', 'rv']
Run_Date_List = ["Dec142013", "Nov302013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]

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
Ax.set_xlabel(r"$\overline{w^{\prime}\theta^{\prime}}$", fontsize=20)
#Ax.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
Ax.set_ylabel(r"$z/z_{g}$", fontsize=20)

plt.xlim(-0.2, 1)
#Ax.set_xticks([-.01, 0, .2])
#plt.ylim(100, 1500)
plt.ylim(0, 1.2)


Ax1 = Fig1.add_subplot(132)
#Ax1.set_title( r'$Scaled \ \frac{\partial \theta}{\partial z}$', fontsize=20)
Ax1.set_xlabel( r'$+$', fontsize=20)
#Ax1.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
#Ax1.set_xlabel(r"$\frac{\partial \theta}{\partial z}$", fontsize=20)
#Ax1.set_ylabel(r"$z$", fontsize=20)
#Ax1.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#start, end = -.025, .025
#start, end = -1, 2.5
#Ax1.set_xticks([0, .005])
#Ax1.set_ylabel(r"$z$", fontsize=20)
#plt.xlim(-.025, .025)
plt.xlim(0, 2)
#plt.ylim(100, 1500)
plt.ylim(0, 1.2)

Ax2 = Fig1.add_subplot(133)
#Ax2.set_title(r"$\overline{w^{'} \theta^{'}}$", fontsize=20)
#Ax2.set_title(r"$Scaled \ \overline{w^{'} \theta^{'}}$", fontsize=20)
Ax2.set_xlabel(r"$-$", fontsize=20)
#Ax2.set_xlabel(r"$\frac{\overline{w^{'}\theta^{'}}}{\overline{w^{'}\theta^{'}}_{0}}$", fontsize=20)
#start, end = -.06, .14
#start, end = -.4, 1.2
#Ax2.set_xticks([0])

#Ax2.set_ylabel(r"$z$", fontsize=20)
#Ax2.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#plt.ylim(100, 1500)
#plt.xlim(-.006, .004)
plt.xlim(-1,0)
plt.ylim(0,1.2)
dump_time_list, Times = Make_Timelists(1, 1800, 28800)
#dump_time_list, Times = dump_time_list[12:], Times[12:]  
theta_file_list = ["/tera/users/nchaparr/"+date+"/data/theta_bar"+ dump_time for dump_time in dump_time_list]
press_file_list = ["/tera/users/nchaparr/"+date+"/data/press"+ dump_time for dump_time in dump_time_list]
flux_file_list = ["/tera/phil/nchaparr/python/Plotting/"+date+"/data/flux_quads_test"+ dump_time for dump_time in dump_time_list]
flux_file_list1 = ["/tera/users/nchaparr/"+date+"/data/wvelthetapert"+ dump_time for dump_time in dump_time_list]
height_file = "/tera/users/nchaparr/"+date+"/data/heights0000018000"
AvProfVars = np.genfromtxt("/tera/users/nchaparr/"+date+"/data/AvProfLims")

#loop over text files files
for i in range(len(theta_file_list)):    
    theta = np.genfromtxt(theta_file_list[i])
    height = np.genfromtxt(height_file)
    press = np.genfromtxt(press_file_list[i])
    rhow = nc.calc_rhow(press, height, theta[0])
    print(rhow[0])
    flux_quads = np.genfromtxt(flux_file_list[i])
    #flux = np.genfromtxt(flux_file_list1[i])
    upwarm=flux_quads[:,0]
    upwarm[0] = np.nan
    downwarm=flux_quads[:,1]
    downwarm[0]=np.nan
    upcold=flux_quads[:,2]
    upcold[0] = np.nan
    downcold=flux_quads[:,3]
    downcold[0] = np.nan
    flux=flux_quads[:,4]
    flux[0] = np.nan
    #Now for the gradients
    dheight = np.diff(height)
    #print dheight
    dtheta = np.diff(theta)      
    dthetadz = np.divide(dtheta, dheight)
           
    element0 = np.array([0])
    dthetadz=np.hstack((dthetadz, element0))/gamma
        
    #only need up to 2500meters
    top_index = np.where(abs(2545 - height) < 40.)[0][0]
    
    if date=="Nov302013":
        h_time_index=(i+1)*2-1
    else:
        h_time_index=(i+1)*3-1
    
    #where gradient is max, and flux is min
    #print AvProfVars[:,1].shape, height.shape
    scaled_height = [1.0*h/AvProfVars[h_time_index,1] for h in height]

    upwarm = np.multiply(upwarm, rhow)*1004.0/sfc_flx
    downwarm = np.multiply(downwarm, rhow)*1004.0/sfc_flx
    upcold = np.multiply(upcold, rhow)*1004.0/sfc_flx
    downcold = np.multiply(downcold, rhow)*1004.0/sfc_flx
    flux = np.multiply(flux, rhow)*1004.0/sfc_flx
    

    if np.mod(i+1, 1) == 0:
    #if i > 45 and i < 49:
        
        upwarm[0] = np.nan
        zeros = np.zeros_like(height)

        Ax.plot(flux, scaled_height, 'k-') #, label = str(Times[i])+'hrs'        
        Ax1.plot(upwarm, scaled_height,'r-')
        Ax1.plot(downcold, scaled_height, 'b-')
        Ax2.plot(upcold, scaled_height,'b-')
        Ax2.plot(downwarm, scaled_height, 'r-')
        
    
array = np.genfromtxt('/tera/phil/nchaparr/python/Pert_Files/snd')
    
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
#Ax1.plot(zeros, scaled_height)#zeros line for reference
#Ax1.plot(gamma, scaled_height)#zeros line for reference
#Ax1.plot(zeros+gamma, height, 'k-')#zeros line for reference
#Ax.plot(zeros, scaled_height)#zeros line for reference
#Ax2.plot(zeros, height)#zeros line for reference 
#plt.legend(loc = 'Lower right', prop={'size':8})
#Ax2.plot(theta_0, scaled_xheight[0:top_index], '--', label = 'Initial Sounding')#"
#plt.xlim(300, 310)
#plt.legend(loc = 'upper right', prop={'size':8})
#print "plotting"
plt.show()
#Fig1.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/theta_flux_profs.png")





    
    
