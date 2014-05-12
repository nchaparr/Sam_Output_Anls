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
        Deltah = points.Deltah_over_h()
        HistVars = points.HistVars()
        AvProfVars = points.AvProfVars()
    #TODO: alternative starting index for Nov302013
        if Run_Date_List[i]=="Jan152014_1":
             Ax3.plot(rinovals[11:29, 1], rinovals[11:29, 9], legend_list[i], label = label_list[i])
        else:
             #print len(Times[11:]), rinovals[11:, 7].shape
             Ax3.plot(rinovals[11:, 1], rinovals[11:, 9], legend_list[i], label = label_list[i])
        
#Ax3.plot(np.arange(0, .1, .01)[2:10], .20833*np.arange(0, .1, .01)[2:10], 'k--')
#Ax3.plot(np.arange(0, .1, .01)[2:10], np.arange(0, .1, .01)[2:10]**(3.0/2), 'k--')
#Ax3.plot(Times[11:], Fit, 'b-', label="2nd Order Polyfit")
#Ax3.text(6, 1500, r'$ \frac{dh}{dt}  =  %.3f \frac{m}{s} $' %(1.0*M/3600),  fontdict=None, withdash=False, fontsize = 15)
#Ax3.set_ylim(0, 2500)
Ax3.legend(loc = 'upper right', prop={'size': 10}, numpoints=1)
#Ax3.set_title(r'$\Delta h (Flux)\ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$Scaled \ Time \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\frac{\Delta h}{h} \ vs \ Ri^{-1}$', fontsize=20)
Ax3.set_xlabel(r'$Ri^{-1}$', fontsize=20)
#Ax3.set_title(r'$\Delta \theta \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\overline{\theta} \ vs \ Time$', fontsize=20)
#Ax3.set_xlabel(r"$\frac{Time}{\tau}$", fontsize=20)
Ax3.set_ylabel(r"$pi3_{Douw}$", fontsize=20)
#Ax3.set_ylabel(r"$\frac{\Delta h}{h}$", fontsize=20)
#Ax3.set_ylabel(r"$\frac{w_{e}}{w^{*}}$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta h (m)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta \theta (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\overline{ \theta} (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta h \ (m)$", fontsize=20)
#Ax3.set_ylabel(r"$z \ (m)$", fontsize=20)
#Ax3.set_xlabel(r"$Time \ (hrs)$", fontsize=20)
#Ax3.set_xlabel(r"$\gamma \frac{\Delta h}{\Delta \theta}$", fontsize=20)
#Ax3.set_ylabel(r"$h \ (m)$", fontsize=20)
#plt.ylim(0, 1.4)
plt.show()





    
    
