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
import pandas as pd

"""

   For plotting the (scaled) temperature gradient and flux profiles.

"""

Fig1 = plt.figure(1)
Fig1.clf()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


Ax = Fig1.add_subplot(111)
#Ax.set_title( r'$\theta$', fontsize=20)
#Ax1.set_title( r'$\frac{\partial \theta}{\partial z}$', fontsize=20)
#Ax1.set_xlabel(r"$\frac{\frac{\partial \theta}{\partial z}}{\gamma}$", fontsize=20)
#Ax.set_xlabel(r"$\frac{ \partial \overline{\theta} }{\partial z} / \gamma$", fontsize=25)
Ax.set_xlabel(r"$\frac{\partial \overline{\theta}}{\partial z}$ Kkm$^{-1}$", fontsize=35)
#Ax.set_xlabel(r"$\overline{w^{,}\theta^{,}}$", fontsize=20)
#Ax1.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
Ax.set_ylabel(r"$\frac{z}{h}$", fontsize=40)
#plt.xlim(0, 1)
plt.xlim(-.0002, .02)
#plt.ylim(50, 950)
plt.ylim(0.2, 1.4)


dump_time_list, Times = Make_Timelists(1, 600, 28800)
marker_list=['ko', 'kv', 'yo', 'y*', 'ro', 'yv', 'rv']
legend_list=["100/5", "100/10", "60/5", "60/2.5", "150/5", "60/10", "150/10"]
run_name_list = ["Nov302013", "Dec142013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]
flux_list = [100, 100, 60, 60, 150, 60, 150]
gamma_list=[.005, .01, .005, .0025, .005, .01, .01]
#choose a dumptime
#dump_time, Time = dump_time_list[time_index], Times[time_index]
##print Time
dump_time = "0000010800"
dump_time_index=29
dump_time_index0=19


theta_file_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+run_name+"/data/theta_bar"+ dump_time for run_name in run_name_list]
press_file_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+run_name+"/data/press"+ dump_time for run_name in run_name_list]
flux_file_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+run_name+"/data/wvelthetapert"+ dump_time for run_name in run_name_list]
height_file_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+run_name+"/data/heights0000000600" for run_name in run_name_list]
AvProfVars_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+run_name+"/data/AvProfLims" for run_name in run_name_list]
invrinos_list = ["/newtera/tera/phil/nchaparr/python/Plotting/"+run_name+"/data/invrinos" for run_name in run_name_list]


with pd.HDFStore('paper_table.h5','r') as store:
     print(store.keys())
     df_cases=store.get('cases')

group_attributes={}
time_full=np.linspace(600,28800,48)
time_nov302013=np.linspace(900,28800,32)
history='written 2015/8/5 by plot_theta_profs.py  9542a821e'        
group_attributes['/']=dict(time600=time_full,time900=time_nov302013,history=history)

with pd.HDFStore('all_profiles.h5','w') as store:
    #loop over text files files
    store.put('df_overview',df_cases,format='table')
    for i in range(len(theta_file_list)):
        run_name = run_name_list[i]
        theta = np.genfromtxt(theta_file_list[i])
        print('theta shape: ',theta.shape)
        height = np.genfromtxt(height_file_list[i])
        df_prof=pd.DataFrame(height,columns=['height'])
        df_prof['theta']=theta
        gamma = gamma_list[i]

        press = np.genfromtxt(press_file_list[i])
        rhow = nc.calc_rhow(press, height, theta[0])
        df_prof['press']=press
        df_prof['rho']=rhow
        wvelthetapert = np.genfromtxt(flux_file_list[i])
        wvelthetapert[0] = np.nan
        df_prof['wvelthetapert']=wvelthetapert
        AvProfVars = np.genfromtxt(AvProfVars_list[i])
        invrinosVars = np.genfromtxt(invrinos_list[i])
        print(run_name,' invrinos ',invrinosVars.shape)
        if AvProfVars.shape[0]==32:
            big_array=np.empty([48,8],dtype=np.float64)
            big_array[:,:]=np.nan
            big_array[0:32,:]=AvProfVars[:,:]
            AvProfVars=big_array[:,:]
        if invrinosVars.shape[0]==32:
            big_array=np.empty([48,10],dtype=np.float64)
            big_array[:,:]=np.nan
            big_array[0:32,:]=invrinosVars[:,:]
            invrinosVars=big_array[:,:]
        print('case: ',run_name,' shape: ',AvProfVars.shape)
        columns=['h0','h','h1','zf0','zf','zf1','deltatheta','mltheta']
        df_lims=pd.DataFrame(AvProfVars,columns=columns)
        columns=['rino', 'invrino', 'wstar', 'S', 'tau', 'mltheta', 'deltatheta', 
                 'pi3', 'pi4', 'thetastar']
        df_rinos=pd.DataFrame(invrinosVars,columns=columns)
       #Now for the gradients
        dheight = np.diff(height)
        dtheta = np.diff(theta)      
        dthetadz = np.divide(dtheta, dheight)        
        element0 = np.array([0])
        dthetadz=np.hstack((element0, 1.0*dthetadz)) #*1.0/gamma
        df_prof['dthetadz']=dthetadz
        node_name='cases/{}/df_lims'.format(run_name)
        store.put(node_name,df_lims,format='table')
        node_name='cases/{}/df_prof'.format(run_name)
        store.put(node_name,df_prof,format='table')
        node_name='cases/{}/df_rinos'.format(run_name)
        store.put(node_name,df_rinos,format='table')


        #only need up to 2500meters
        top_index = np.where(abs(1670 - height) < 40.)[0][0]

        #where gradient is max, and flux is min
        ##print AvProfVars[:,1].shape, height.shape

        if run_name == "Nov302013":
            h1 = AvProfVars[dump_time_index0, 1]
        else:
            h1 = AvProfVars[dump_time_index, 1]

        h_index=np.where(dthetadz - np.amax(dthetadz[:top_index])==0)[0][0]
        h=height[h_index]    
        scaled_height = [1.0*ht/h for ht in height]

        #print h1, h_index, height[h_index]

        fluxes = np.multiply(wvelthetapert, rhow)*1004.0/flux_list[i]

        Ax.plot(dthetadz, scaled_height, marker_list[i], label = legend_list[i], markersize=10) #, 

zeros = np.zeros_like(height)
Ax.plot(zeros+.01, scaled_height, 'k-')
Ax.plot(zeros+.005, scaled_height, 'k-')
Ax.plot(zeros+.005, scaled_height, 'k-')
Ax.plot(zeros+.0025, scaled_height, 'k-')
Ax.plot(zeros+.0002, scaled_height, 'k-')
#Ax.legend(numpoints=1, loc = 'lower right', prop={'size':14})
Ax.set_xticks([0.0002, .0025, .005, .010])
Ax.set_xticklabels([0.2, 2.5, 5, 10])
Ax.tick_params(axis="both", labelsize=25)
plt.tight_layout()
plt.show()
#Fig1.savefig('/tera/phil/nchaparr/python/Plotting/Dec252013/pngs/theta_profs2hrs.png')
#Fig1.savefig('/tera/phil/nchaparr/python/Plotting/Dec252013/pngs/flux_profs2hrs.png')





    
    
