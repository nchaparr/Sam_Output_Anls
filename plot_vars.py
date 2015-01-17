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
   Plots 

"""

Fig1 = plt.figure(1)
Fig1.clf()

Ax1 = Fig1.add_subplot(111)
#Ax1.set_title(r'Scaled Root Mean Velocity Squared - 2 hours')
Ax1.set_xlabel(r"$\frac{\sqrt[]{u^{,2}}}{w^{*}}$", fontsize=20)
Ax1.set_ylabel(r'$\frac{z}{h}$', fontsize=20)

#plt.xlim(300, 310)

dump_time_list, Times = Make_Timelists(1, 3600, 28800)
dump_time1_list = ['0000003600', '0000007200', '0000010800', '0000014400', '0000018000', '0000021600', '0000025200', '0000028800']
flux_file_list = ["/tera/phil/nchaparr/python/Plotting/Mar52014/data/wvelthetapert"+ dump_time for dump_time in dump_time_list]
u_file_list = ["/tera/phil/nchaparr/python/Plotting/Mar52014/data/rootmeanusq"+ dump_time for dump_time in dump_time_list] 
v_file_list = ["/tera/phil/nchaparr/python/Plotting/Mar52014/data/rootmeanvsq"+ dump_time for dump_time in dump_time_list]
w_file_list = ["/tera/phil/nchaparr/python/Plotting/Mar52014/data/rootmeanwsq"+ dump_time for dump_time in dump_time_list]
#flux_quad_file_list = ["/tera/phil/nchaparr/python/Plotting/Sep302013/data/flux_quads"+ dump_time1 for dump_time1 in dump_time1_list]
height_file = "/tera/phil/nchaparr/python/Plotting/Mar52014/data/heights0000028800"

colorlist=['k', 'b', 'c', 'g', 'r', 'm', 'y', '.75']
#loop over text files files
for i in range(len(flux_file_list)):    
     
    wvelthetapert = np.genfromtxt(flux_file_list[i])
    #print flux_file_list[(i+1)*60-1]
    u  = np.genfromtxt(u_file_list[i])
    v  = np.genfromtxt(v_file_list[i])
    w  = np.genfromtxt(w_file_list[i])
    invrinos = np.genfromtxt('/tera/phil/nchaparr/python/Plotting/Mar52014/data/invrinos')
    AvProfLims = np.genfromtxt('/tera/phil/nchaparr/python/Plotting/Mar52014/data/AvProfLims')
    height = np.genfromtxt(height_file)
    zeros = np.zeros_like(height)
    
    if np.mod(i+1, 1) ==0 and i == 2:
        color = colorlist[int(1.0*i/6)]
        wstar = invrinos[(i+1)*6-1,2]
        h = AvProfLims[(i+1)*6-1,2]
        print h
        Ax1.plot(1.0*u/wstar, 1.0*height/h, 'r--', label='u') #line for reference
        Ax1.plot(1.0*v/wstar, 1.0*height/h, 'b*', label='v')
        Ax1.plot(1.0*w/wstar, 1.0*height/h, 'g-', label='w')          
        
        #Ax1.plot(wvelthetapert, height, color, label = str(Times[i])+'hrs')
        #Ax1.plot(zeros, height, color, label = str(Times[i])+'hrs')
     
   
plt.ylim(0,1.5)
plt.xlim(0, 1)
plt.legend(loc = 'upper right', prop={'size':8}, numpoints=1)
plt.show()
Fig1.savefig('/tera/phil/nchaparr/python/Plotting/Mar52014/pngs/rmsvel2.png')






    
    
