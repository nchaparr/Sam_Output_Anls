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

"""plots scaled dhdt vs inverse richardson number based on average profile
   now has option to do least squares, and plot 
"""

dump_time_list, Times = Make_Timelists(1, 600, 28800)
Times = np.array(Times)
HistVars =  np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec252013/data/ELLims_hist")
AvProfVars = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec252013/data/AvProfLims")
rinovals = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec252013/data/invrinos")
#plot the heights vs time

Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)
FitFunc=np.polyfit(Times[11:], AvProfVars[11:, 1], 2, full=False)
Fit = FitFunc[0]*Times[11:]**2 + FitFunc[1]*Times[11:] + FitFunc[2]
dhdt =1.0*(2*FitFunc[0]*Times[11:] + FitFunc[1])/3600
#Fit = FitFunc[0]*Times[120:]**3 + FitFunc[1]*Times[120:]**2 + FitFunc[2]*Times[120:] + FitFunc[3]
#dhdt =1.0*(3*FitFunc[0]*Times[120:]**2 + 2*FitFunc[1]*Times[120:] + FitFunc[2])/3600
deltah = np.subtract(AvProfVars[:, 2], AvProfVars[:, 0])
deltah = np.divide(deltah, AvProfVars[:, 1])
tau = 1.0*rinovals[:,4]/3600
scaled_time = np.divide(Times, tau)
scaled_dhdt = np.divide(dhdt, rinovals[11:, 2])
dhdtinvriplt = np.vstack((rinovals[11:, 1], scaled_dhdt))
dhdtinvriplt = np.transpose(np.vstack((dhdtinvriplt,deltah[11:])))
np.savetxt('/tera/phil/nchaparr/python/Plotting/Dec252013/data/dhdtinvriplt.txt', dhdtinvriplt, delimiter=' ')
Ax3.plot(rinovals[11:, 1], deltah[11:], 'ko')
#Ax3.plot(Times[:], rinovals[:, 1], 'b*', label="")
#Ax3.plot(Times[11:], AvProfVars[11:, 1], 'ko', label = 'Avg Prof')
#Ax3.plot(Times[11:], Fit, 'b-', label="2nd Order Polyfit")
#Ax3.text(6, 1500, r'$ \frac{dh}{dt}  =  %.3f \frac{m}{s} $' %(1.0*M/3600),  fontdict=None, withdash=False, fontsize = 15)
#Ax3.set_ylim(0, 2500)
Ax3.legend(loc = 'upper left', prop={'size': 10})
Ax3.set_title(r'$Scaled \ \Delta h \ vs\ Ri^{-1}$', fontsize=20)
Ax3.set_ylabel(r"$\frac{\Delta h}{h}$", fontsize=20)
Ax3.set_xlabel(r"$Ri^{-1}$", fontsize=20)
#plt.ylim(.015, .04)
plt.show()





    
    
