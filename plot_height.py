import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
from nchap_class import *
import nchap_fun as nc
from matplotlib.lines import Line2D

from matplotlib import rcParams
rcParams.update({'font.size': 10})

"""
   plots heights, rinos, etc
  
"""

dump_time_list, Times = Make_Timelists(1, 600, 28800)
Times = np.array(Times)
dump_time_list0, Times0 = Make_Timelists(1, 900, 28800)
Times0=np.array(Times0)
#plot the heights vs time
#print Line2D.markers
Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)

#This needs to be included, for scaled time vs height plot
#tau = 1.0*rinovals5[:,4]/3600
#scaled_time = np.divide(Times, tau)


#Main Part -- pulling points and plotting them
label_list = ['100/10', '100/5', '60/5', '60/2.5', '150/5', '60/10', '150/10']
legend_list = ['kv', 'ko', 'yo', 'y*', 'ro', 'yv', 'rv']
Run_Date_List = ["Dec142013", "Nov302013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]


for i in range(len(label_list)):
    if i<99:
        points = For_Plots(Run_Date_List[i])
        rinovals = points.rinovals()
        Deltah = points.Deltah_over_h()
        HistVars = points.HistVars()
        AvProfVars = points.AvProfVars()
        #dhdtinvriplt=For_Plots.get_dhdt(Times, 11, 47)  TODO: test this method for the scaled we vs invri plot
        if Run_Date_List[i]=="Jan152014_1":
              #print len(Times[11:]), AvProfVars[11:, 0].shape
              #Ax3.plot(Times[11:],  np.divide(AvProfVars[11:, 2], AvProfVars[11:, 1]), legend_list[i], label = label_list[i])
             rinovals[16:21, 1] = np.nan
             Deltah[16:21] = np.nan
             Ax3.plot(Times[11:29], np.divide(AvProfVars[11:29, 0], AvProfVars[11:29, 1]), legend_list[i])
             Ax3.plot(Times[11:29], np.divide(AvProfVars[11:29, 2], AvProfVars[11:29, 1]), legend_list[i], label = label_list[i])
        elif Run_Date_List[i]=="Nov302013":
             #
             #Ax3.plot(Times[11:],  np.divide(AvProfVars[11:, 0], AvProfVars[11:, 1]), legend_list[i], label = label_list[i])
             Ax3.plot(Times0[7:], np.divide(AvProfVars[7:, 0], AvProfVars[7:, 1]), legend_list[i])
             Ax3.plot(Times0[7:], np.divide(AvProfVars[7:, 2], AvProfVars[7:, 1]), legend_list[i], label = label_list[i])
        elif Run_Date_List[i] =="Mar12014":
             rinovals[12:17, 1] = np.nan
             Deltah[12:17] = np.nan       
             Ax3.plot(Times[11:], np.divide(AvProfVars[11:, 0], AvProfVars[11:, 1]), legend_list[i])
             Ax3.plot(Times[11:], np.divide(AvProfVars[11:, 2], AvProfVars[11:, 1]), legend_list[i], label = label_list[i])
        else:
             Ax3.plot(Times[11:], np.divide(AvProfVars[11:, 0], AvProfVars[11:, 1]), legend_list[i])
             Ax3.plot(Times[11:], np.divide(AvProfVars[11:, 2], AvProfVars[11:, 1]), legend_list[i], label = label_list[i])

#Ax3.plot(np.arange(0, .1, .01)[2:10], .20833*np.arange(0, .1, .01)[2:10], 'k--')
#Ax3.plot(np.arange(0, .1, .01)[2:10], np.arange(0, .1, .01)[2:10]**(3.0/2), 'k--')
#Ax3.plot(Times[11:], Fit, 'b-', label="2nd Order Polyfit")
#Ax3.text(6, 1500, r'$ \frac{dh}{dt}  =  %.3f \frac{m}{s} $' %(1.0*M/3600),  fontdict=None, withdash=False, fontsize = 15)
#Ax3.set_ylim(0, 2500)
Ax3.legend(loc = 'lower left', prop={'size': 10}, numpoints=1)
#Ax3.set_title(r'$\Delta h (Flux)\ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$Scaled \ Time \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\frac{\Delta h}{h} \ vs \ Ri^{-1}$', fontsize=20)
Ax3.set_xlabel(r'$Time \ (hrs)$', fontsize=20)
#Ax3.set_title(r'$\Delta \theta \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\overline{\theta} \ vs \ Time$', fontsize=20)
#Ax3.set_xlabel(r"$\frac{Time}{\tau}$", fontsize=20)
#Ax3.set_ylabel(r"$pi3_{Douw}$", fontsize=20)
Ax3.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#Ax3.set_ylabel(r"$\frac{w_{e}}{w^{*}}$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta h (m)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta \theta (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\overline{ \theta} (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta h \ (m)$", fontsize=20)
#Ax3.set_ylabel(r"$z \ (m)$", fontsize=20)
#Ax3.set_xlabel(r"$Time \ (hrs)$", fontsize=20)
#Ax3.set_xlabel(r"$\gamma \frac{\Delta h}{\Delta \theta}$", fontsize=20)
#Ax3.set_ylabel(r"$h \ (m)$", fontsize=20)
plt.ylim(0, 1.5)
#plt.xlim(.02, .1)
plt.show()





    
    
