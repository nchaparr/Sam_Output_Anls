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

   For plotting the height definitions.

"""

date = "Mar52014"
sfc_flx = 150
gamma = .01

Fig1 = plt.figure(1)
Fig1.clf()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


Ax = Fig1.add_subplot(131)
#Ax.set_title( r'$\theta$', fontsize=20)
#Ax.set_title( r'$\frac{\partial \theta}{\partial z}$', fontsize=20)
#Ax.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
Ax.set_xlabel(r"$\overline{\theta} \ (K)$", fontsize=20)
#Ax.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
Ax.set_ylabel(r"$z \ (m)$", fontsize=20)
plt.xlim(300, 320)
plt.ylim(100, 1500)
#plt.ylim(0.1, 1.4)


Ax1 = Fig1.add_subplot(132)
#Ax1.set_title( r'$Scaled \ \frac{\partial \theta}{\partial z}$', fontsize=20)
#Ax1.set_title( r'$\frac{\partial \theta}{\partial z}$', fontsize=20)
#Ax1.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
Ax1.set_xlabel(r"$\frac{\partial \theta}{\partial z} \ (\frac{K}{m})$", fontsize=20)
#Ax1.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
start, end = -.025, .025
#start, end = -1, 2.5
Ax1.set_xticks(np.arange(start, end, 1.0*(end-start)/4))
#Ax1.set_ylabel(r"$z$", fontsize=20)
plt.xlim(-.025, .025)
#plt.xlim(-1, 2.5)
plt.ylim(100, 1500)
#plt.ylim(0.1, 1.4)

Ax2 = Fig1.add_subplot(133)
#Ax2.set_title(r"$\overline{w^{'} \theta^{'}}$", fontsize=20)
#Ax2.set_title(r"$Scaled \ \overline{w^{'} \theta^{'}}$", fontsize=20)
Ax2.set_xlabel(r"$\overline{w^{'}\theta^{'}} \ (\frac{m}{s}K)$", fontsize=20)
#Ax2.set_xlabel(r"$\frac{\overline{w^{'}\theta^{'}}}{\overline{w^{'}\theta^{'}}_{0}}$", fontsize=20)
start, end = -.08, .14
plt.xlim(-.08, .14)
#start, end = -.6, 1.2
Ax2.set_xticks(np.arange(start, end, 1.0*(end-start)/5))

#Ax2.set_ylabel(r"$z$", fontsize=20)
#Ax2.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
plt.ylim(100, 1500)
#plt.xlim(-.06, .14)
#plt.xlim(-.4, 1.2)
#plt.ylim(0.1, 1.4)
dump_time_list, Times = Make_Timelists(1, 600, 28800)
 
theta_file_list = ["/tera/phil/nchaparr/python/Plotting/"+date+"/data/theta_bar"+ dump_time for dump_time in dump_time_list]
press_file_list = ["/tera/phil/nchaparr/python/Plotting/"+date+"/data/press"+ dump_time for dump_time in dump_time_list]
flux_file_list = ["/tera/phil/nchaparr/python/Plotting/"+date+"/data/wvelthetapert"+ dump_time for dump_time in dump_time_list]
height_file = "/tera/phil/nchaparr/python/Plotting/"+date+"/data/heights0000000600"
AvProfVars = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/AvProfLims")

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

    
    #if np.mod(i+1, 6) == 0:
    if i == 27:
        h0=AvProfVars[i,0]
        h=AvProfVars[i,1]
        h1=AvProfVars[i,2]

        z_f0=AvProfVars[i,3]
        z_f=AvProfVars[i,4]
        z_f1=AvProfVars[i,5]

        h0_index=np.where(height==h0)[0]
        h_index=np.where(height==h)[0]
        h1_index=np.where(height==h1)[0]
    
        z_f0_index=np.where(height==z_f0)[0]
        z_f_index=np.where(height==z_f)[0]
        z_f1_index=np.where(height==z_f1)[0]
    
        fluxes[0] = np.nan
        zeros = np.zeros_like(height)

        Ax.plot(theta, height, 'k-', label = r"$\overline{\theta}$") # 
        #Ax.plot([theta[], theta[]], [height[], height[]], '')

        Ax1.plot(1.0*dthetadz, height, 'k-') #, label = str(Times[i])+'hrs'
        Ax1.plot([dthetadz[h0_index]-.02, dthetadz[h0_index]+.023], [h0, h0], 'k--')
        Ax1.plot([dthetadz[h0_index]-.02, dthetadz[h0_index]+.023], [h, h], 'k-')
        Ax1.plot([dthetadz[h0_index]-.02, dthetadz[h0_index]+.023], [h1, h1], 'k--')

        Ax1.text(dthetadz[h0_index]-.02, h1, r"$h_{1}$", size=20)
        Ax1.text(dthetadz[h0_index]-.02, h, r"$h$", size=20)
        Ax1.text(dthetadz[h0_index]-.02, h0, r"$h_{0}$", size=20)
                
        
        Ax2.plot(wvelthetapert, height, 'k-') #, label = str(Times[i])+'hrs'    
        Ax2.plot([wvelthetapert[z_f0_index]-.05, wvelthetapert[z_f0_index]+.1], [z_f0, z_f0], 'k--')
        Ax2.plot([wvelthetapert[z_f0_index]-.05, wvelthetapert[z_f0_index]+.1], [z_f1, z_f1], 'k--')
        Ax2.plot([wvelthetapert[z_f0_index]-.05, wvelthetapert[z_f0_index]+.1], [z_f, z_f], 'k-')

        Ax2.text(dthetadz[z_f0_index]+.10, z_f0, r"$z_{f0}$", size=20)
        Ax2.text(dthetadz[z_f0_index]+.10, z_f1, r"$z_{f1}$", size=20)
        Ax2.text(dthetadz[z_f0_index]+.10, z_f, r"$z_{f}$", size=20)
       
        
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

Ax.plot(theta_0 -.2, height[0:top_index], 'k--', label = r"$\overline{\theta}_{0}$") #, 
#Ax.plot([theta_0[], theta[]], [height[], height[]], '--')

#theAx.text(300, 1500, '',  fontdict=None, withdash=False)
#theAx.text(300, 1400, '',  fontdict=None, withdash=False)

Ax1.plot(zeros, height, 'k-')#zeros line for reference
#Ax1.plot(gamma, height)#zeros line for reference
Ax1.plot(zeros+gamma, height, 'k-')#zeros line for reference
#Ax2.plot(zeros, height)#zeros line for reference
Ax2.plot(zeros, height, 'k-')#zeros line for reference 
#Ax2.plot(theta_0, scaled_xheight[0:top_index], '--', label = 'Initial Sounding')#"
#plt.xlim(300, 310)
Ax.legend(loc = 'upper left', prop={'size':14})
plt.show()
#Fig1.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/theta_flux_profs.png")





    
    
