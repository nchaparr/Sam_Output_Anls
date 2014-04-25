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
FitFunc=np.polyfit(Times[11:], AvProfVars5[11:, 1], 2, full=False)
Fit = FitFunc[0]*Times[11:]**2 + FitFunc[1]*Times[11:] + FitFunc[2]
dhdt =1.0*(2*FitFunc[0]*Times[11:] + FitFunc[1])/3600

#Fit = FitFunc[0]*Times[120:]**3 + FitFunc[1]*Times[120:]**2 + FitFunc[2]*Times[120:] + FitFunc[3]
#dhdt =1.0*(3*FitFunc[0]*Times[120:]**2 + 2*FitFunc[1]*Times[120:] + FitFunc[2])/3600

#Not sure I need this, doing it already above
deltah = np.subtract(AvProfVars5[:, 2], AvProfVars5[:, 0])
deltah = np.divide(deltah, AvProfVars5[:, 1])
tau = 1.0*rinovals5[:,4]/3600
scaled_time = np.divide(Times, tau)

#saving the scaled we vs invri plot points
scaled_dhdt = np.divide(dhdt, rinovals5[11:, 2])
dhdtinvriplt = np.vstack((rinovals5[11:, 1], scaled_dhdt))
dhdtinvriplt = np.transpose(np.vstack((dhdtinvriplt,deltah[11:])))
np.savetxt('/tera/phil/nchaparr/python/Plotting/Mar52014/data/dhdtinvriplt.txt', dhdtinvriplt, delimiter=' ')

#Main Part -- pulling points and plotting them
label_list = ['100/10', '100/5', '60/5', '60/2.5', '150/5', '60/10', '150/10']
legend_list = ['kv', 'ko', 'yo', 'y*', 'ro', 'yv', 'rv']
Run_Date_List = ["Dec142013", "Nov302013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]


for i in range(len(label_list)):
    points = For_Plots(Run_Date_List[i])
    rinovals = points.rinovals()
    Deltah = points.Deltah_over_h()
    HistVars = points.HistVars()
    AvProfVars = points.AvProfVars()
    #TODO: alternative starting index for Nov302013
    Ax3.plot(rinovals[11:,8], Deltah[11:], legend_list[i], label = label_list[i])

#Ax3.plot(rinovals0[7:,8], Deltah00[7:], 'ko', label = '100/5')
#Ax3.plot(rinovals1[11:,8], Deltah01[11:], 'yo', label = '60/5')
#Ax3.plot(rinovals2[11:,8], Deltah02[11:], 'y*', label = '60/2.5')
#Ax3.plot(rinovals3[11:29,8], Deltah03[11:29], 'ro', label = '150/5')
#Ax3.plot(rinovals4[11:,8], Deltah04[11:], 'yv', label = '60/10')
#Ax3.plot(rinovals5[11:,8], Deltah05[11:], 'rv', label = '150/10')

#Ax3.plot(Times[11:], HistVars[11:, 1], 'b*', label="h Dist")
#Ax3.plot(rinovals[:,1], deltah[:], 'ko', label = '100/10')
#print len(Times0[7:]), Deltah00[:].shape
#Ax3.plot(Times[11:], np.divide(AvProfVars[11:,0], AvProfVars[11:,1]), linestyle='none', marker = 'v', color = 'k', label = '100/10')
#Ax3.plot(Times[11:], np.divide(AvProfVars[11:,2], AvProfVars[11:,1]), linestyle='none', marker = 'v', color = 'k')
#Ax3.plot(Times0[7:], np.divide(AvProfVars0[7:,0], AvProfVars0[7:,1]), linestyle='none', marker = 'o', color = 'k', label = '100/5')
#Ax3.plot(Times0[7:], np.divide(AvProfVars0[7:,2], AvProfVars0[7:,1]), linestyle='none', marker = 'o', color = 'k')
#Ax3.plot(Times[11:], np.divide(AvProfVars1[11:,0], AvProfVars1[11:,1]), linestyle='none', marker = 'o', color = 'y', label = '60/5')
#Ax3.plot(Times[11:], np.divide(AvProfVars1[11:,2], AvProfVars1[11:,1]), linestyle='none', marker = 'o', color = 'y')
#Ax3.plot(Times[11:], np.divide(AvProfVars2[11:,0], AvProfVars2[11:,1]), linestyle='none', marker = '*', color = 'y', label = '60/2.5')
#Ax3.plot(Times[11:], np.divide(AvProfVars2[11:,2], AvProfVars2[11:,1]), linestyle='none', marker = '*', color = 'y')
#Ax3.plot(Times[11:29], np.divide(AvProfVars3[11:29,0], AvProfVars3[11:29,1]), linestyle='none', marker = 'o', color = 'r', label = '150/5')
#Ax3.plot(Times[11:29], np.divide(AvProfVars3[11:29,2], AvProfVars3[11:29,1]), linestyle='none', marker = 'o', color = 'r')
#Ax3.plot(Times[11:], np.divide(AvProfVars4[11:,0], AvProfVars4[11:,1]), linestyle='none', marker = 'v', color = 'y', label = '60/10')
#Ax3.plot(Times[11:], np.divide(AvProfVars4[11:,2], AvProfVars4[11:,1]), linestyle='none', marker = 'v', color = 'y')
#Ax3.plot(Times[11:], np.divide(AvProfVars5[11:,0], AvProfVars5[11:,1]), linestyle='none', marker = 'v', color = 'r', label = '150/10')
#Ax3.plot(Times[11:], np.divide(AvProfVars5[11:,2], AvProfVars5[11:,1]), linestyle='none', marker = 'v', color = 'r')
#Ax3.plot(np.arange(0, .1, .01)[2:10], .20833*np.arange(0, .1, .01)[2:10], 'k--')
#Ax3.plot(np.arange(0, .1, .01)[2:10], np.arange(0, .1, .01)[2:10]**(3.0/2), 'k--')
#Ax3.plot(Times[11:], Fit, 'b-', label="2nd Order Polyfit")
#Ax3.text(6, 1500, r'$ \frac{dh}{dt}  =  %.3f \frac{m}{s} $' %(1.0*M/3600),  fontdict=None, withdash=False, fontsize = 15)
#Ax3.set_ylim(0, 2500)
Ax3.legend(loc = 'lower right', prop={'size': 10}, numpoints=1)
#Ax3.set_title(r'$\Delta h (Flux)\ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$Scaled \ Time \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\frac{\Delta h}{h} \ vs \ Ri^{-1}$', fontsize=20)
#Ax3.set_title(r'$Ri^{-1} \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\Delta \theta \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\overline{\theta} \ vs \ Time$', fontsize=20)
#Ax3.set_xlabel(r"$\frac{Time}{\tau}$", fontsize=20)
#Ax3.set_ylabel(r"$\frac{\Delta h}{h}$", fontsize=20)
Ax3.set_ylabel(r"$\frac{\Delta h}{h}$", fontsize=20)
#Ax3.set_ylabel(r"$\frac{w_{e}}{w^{*}}$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta h (m)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta \theta (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\overline{ \theta} (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta h \ (m)$", fontsize=20)
#Ax3.set_ylabel(r"$z \ (m)$", fontsize=20)
#Ax3.set_xlabel(r"$Time (hrs)$", fontsize=20)
Ax3.set_xlabel(r"$\gamma \frac{\Delta h}{\Delta \theta}$", fontsize=20)
#Ax3.set_ylabel(r"$h \ (m)$", fontsize=20)
#plt.ylim(0, 1)
plt.show()





    
    
