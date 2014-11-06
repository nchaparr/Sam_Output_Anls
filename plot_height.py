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
Times0 = np.array(Times0)
#plot the heights vs time
#print Line2D.markers
Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)

#Getting w_{e} from a polyfit to the height vs time plot
#FitFunc=np.polyfit(Times[11:], AvProfVars5[11:, 1], 2, full=False)
#Fit = FitFunc[0]*Times[11:]**2 + FitFunc[1]*Times[11:] + FitFunc[2]
#dhdt =1.0*(2*FitFunc[0]*Times[11:] + FitFunc[1])/3600

#Fit = FitFunc[0]*Times[120:]**3 + FitFunc[1]*Times[120:]**2 + FitFunc[2]*Times[120:] + FitFunc[3]
#dhdt =1.0*(3*FitFunc[0]*Times[120:]**2 + 2*FitFunc[1]*Times[120:] + FitFunc[2])/3600

#Not sure I need this, doing it already above
#deltah = np.subtract(AvProfVars5[:, 2], AvProfVars5[:, 0])
#deltah = np.divide(deltah, AvProfVars5[:, 1])
#tau = 1.0*rinovals5[:,4]/3600
#scaled_time = np.divide(Times, tau)

#This is an important step -- perhaps should be a function?
#saving the scaled we vs invri plot points
#scaled_dhdt = np.divide(dhdt, rinovals5[11:, 2])
#dhdtinvriplt = np.vstack((rinovals5[11:, 1], scaled_dhdt))
#dhdtinvriplt = np.transpose(np.vstack((dhdtinvriplt,deltah[11:])))
#np.savetxt('/tera/phil/nchaparr/python/Plotting/Mar52014/data/dhdtinvriplt.txt', dhdtinvriplt, delimiter=' ')

#Main Part -- pulling points and plotting them
label_list = ['100/10', '100/5', '60/5', '60/2.5', '150/5', '60/10', '150/10']
legend_list = ['kv', 'ko', 'yo', 'y*', 'ro', 'yv', 'rv']
Run_Date_List = ["Dec142013", "Nov302013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]


for i in range(len(label_list)):
    if i<99:
        points = For_Plots(Run_Date_List[i])
        rinovals = points.rinovals()
        Deltah = points.Deltahf_over_zf()
        HistVars = points.HistVars()
        AvProfVars = points.AvProfVars()
        if Run_Date_List[i] == "Nov302013":
             #Deltah[13] = np.nan
             #Deltah[15:17] = np.nan
             #Deltah[24:26] = np.nan
             points.Get_and_save_dhdt(Times0[7:], AvProfVars[7:, 4], rinovals[7:, 2], rinovals[7:, 1])
             scaled_we_plot = points.scaled_we_plot()
             Ax3.loglog(scaled_we_plot[0, :], scaled_we_plot[1, :], legend_list[i], label = label_list[i], markersize=12)
             #Ax3.loglog(rinovals[7:, 1], Deltah[7:], legend_list[i], label = label_list[i])
        elif Run_Date_List[i] == "Jan152014_1":
    #TODO: alternative starting index for Nov302013
             #Deltah[16:21] = np.nan
             #print Deltah
             points.Get_and_save_dhdt(Times[11:29], AvProfVars[11:29, 4], rinovals[11:29, 2], rinovals[11:29, 1])
             scaled_we_plot = points.scaled_we_plot()
             Ax3.loglog(scaled_we_plot[0, :], scaled_we_plot[1, :], legend_list[i], label = label_list[i], markersize=12)            
             #Ax3.loglog(rinovals[11:29, 1], Deltah[11:29], legend_list[i], label = label_list[i])
        elif Run_Date_List[i] == "Mar12014":
    #TODO: alternative starting index for Nov302013
             #Deltah[11:17] = np.nan
             #print Deltah
             points.Get_and_save_dhdt(Times[11:], AvProfVars[11:, 4], rinovals[11:, 2], rinovals[11:, 1])
             scaled_we_plot = points.scaled_we_plot()
             Ax3.loglog(scaled_we_plot[0, :], scaled_we_plot[1, :], legend_list[i], label = label_list[i], markersize=12)            
                          
             #Ax3.loglog(rinovals[11:29, 1], Deltah[11:29], legend_list[i], label = label_list[i])
        else:
             points.Get_and_save_dhdt(Times[11:], AvProfVars[11:, 4], rinovals[11:, 2], rinovals[11:, 1])
             scaled_we_plot = points.scaled_we_plot()
             Ax3.loglog(scaled_we_plot[0, :], scaled_we_plot[1, :], legend_list[i], label = label_list[i], markersize=12)             
             #Ax3.loglog(rinovals[11:, 1], Deltah[11:], legend_list[i], label = label_list[i])

xes = np.arange(0.033, .053, .0001)
x1es = np.arange(.12, .35, .0001)
ys = 2.2*xes**(1.5)
ys1= .1*x1es**(1)
#Ax3.loglog(xes, ys, 'k--')
Ax3.loglog(x1es, ys1, 'k--')
#Ax3.plot(np.arange(0, .1, .01)[2:10], np.arange(0, .1, .01)[2:10]**(3.0/2), 'k--')
#Ax3.plot(Times[11:], Fit, 'b-', label="2nd Order Polyfit")
Ax3.text(.14, .024, r'$a = -1$',  fontdict=None, withdash=False, fontsize = 30, rotation=18)

#Ax3.text(.035, .024, r'$a = -\frac{3}{2}$',  fontdict=None, withdash=False, fontsize = 30, rotation=21)

#Ax3.set_ylim(0, 2500)
Ax3.legend(loc = 'upper left', prop={'size': 20}, numpoints=1)
#Ax3.set_title(r'$\Delta h (Flux)\ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$Scaled \ Time \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\frac{\Delta h}{h} \ vs \ Ri^{-1}$', fontsize=20)
#Ax3.set_title(r'$Ri^{-1} \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\Delta \theta \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\overline{\theta} \ vs \ Time$', fontsize=20)
#Ax3.set_xlabel(r"$\frac{Time}{\tau}$", fontsize=20)
#Ax3.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#Ax3.set_ylabel(r"$\frac{\Delta h_{f}}{z_{f}}$", fontsize=20)
Ax3.set_ylabel(r"$\frac{w_{e}}{w^{*}}$", fontsize=30)
#Ax3.set_ylabel(r"$\Delta h (m)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta \theta (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\overline{ \theta} (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta h \ (m)$", fontsize=20)
#Ax3.set_ylabel(r"$z \ (m)$", fontsize=20)
Ax3.set_xlabel(r"$Ri_{\delta}^{-1}$", fontsize=30)
Ax3.set_yticks([.2, .4, .6, .8])
Ax3.set_yticklabels([.2, .4, .6, .8])
#Ax3.set_xticks([.08, .12, 1])
Ax3.set_xticklabels([.2, .4, .6, .8, 1])
Ax3.set_xticks([.2, .4, .6, .8, 1])
#Ax3.set_xticklabels([.02, .04, .06, .08, .1])
Ax3.tick_params(axis="both", labelsize=20)
#Ax3.set_xlabel(r"$\gamma \frac{\Delta h}{\Delta \theta}$", fontsize=20)
#Ax3.set_ylabel(r"$h \ (m)$", fontsize=20)
plt.xlim(0.1, 1)
plt.ylim(0.008, .8)
plt.tight_layout()
plt.show()





    
    
