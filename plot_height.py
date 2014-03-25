import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from matplotlib.lines import Line2D

from matplotlib import rcParams
rcParams.update({'font.size': 10})

"""plots scaled dhdt vs inverse richardson number based on average profile
   now has option to do least squares, and plot 
"""

dump_time_list, Times = Make_Timelists(1, 600, 28800)
Times = np.array(Times)
dump_time_list0, Times0 = Make_Timelists(1, 900, 28800)

dhdtplot = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec142013/data/dhdtinvriplt.txt") 
dhdtplot0 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Nov302013/data/dhdtinvriplt.txt")
dhdtplot1 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec202013/data/dhdtinvriplt.txt")
dhdtplot2 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec252013/data/dhdtinvriplt.txt")
dhdtplot3 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Jan152014_1/data/dhdtinvriplt.txt")
dhdtplot4 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Mar12014/data/dhdtinvriplt.txt")
dhdtplot5 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Mar52014/data/dhdtinvriplt.txt")

#HistVars =  np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec142013/data/ELLims_hist")
#HistVars0 =  np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Nov302013/data/ELLims_hist")
#HistVars1 =  np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec202013/data/ELLims_hist")
#HistVars2 =  np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec252013/data/ELLims_hist")
AvProfVars = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec142013/data/AvProfLims")
AvProfVars0 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Nov302013/data/AvProfLims")
AvProfVars1 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec202013/data/AvProfLims")
AvProfVars2 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec252013/data/AvProfLims")
AvProfVars3 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Jan152014_1/data/AvProfLims")
AvProfVars4 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Mar12014/data/AvProfLims")
AvProfVars5 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Mar52014/data/AvProfLims")

rinovals = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec142013/data/invrinos")
rinovals0 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Nov302013/data/invrinos")
rinovals1 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec202013/data/invrinos")
rinovals2 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec252013/data/invrinos")
rinovals3 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Jan152014_1/data/invrinos")
rinovals4 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Mar12014/data/invrinos")
rinovals5 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Mar52014/data/invrinos")

Deltah0 = np.subtract(AvProfVars[:,2], AvProfVars[:,0])
#Deltah0 = np.divide(Deltah0, AvProfVars[:,1])
Deltah00 = np.subtract(AvProfVars0[:,2], AvProfVars0[:,0])
#Deltah00 = np.divide(Deltah00, AvProfVars0[:,1])
Deltah01 = np.subtract(AvProfVars1[:,2], AvProfVars1[:,0])
#Deltah01 = np.divide(Deltah01, AvProfVars1[:,1])
Deltah02 = np.subtract(AvProfVars2[:,2], AvProfVars2[:,0])
#Deltah02 = np.divide(Deltah02, AvProfVars2[:,1])
Deltah03 = np.subtract(AvProfVars3[:,2], AvProfVars3[:,0])
#Deltah03 = np.divide(Deltah03, AvProfVars3[:,1])
Deltah04 = np.subtract(AvProfVars4[:,2], AvProfVars4[:,0])
#Deltah04 = np.divide(Deltah04, AvProfVars4[:,1])
Deltah05 = np.subtract(AvProfVars5[:,2], AvProfVars5[:,0])
#Deltah05 = np.divide(Deltah05, AvProfVars5[:,1])

#plot the heights vs time
#print Line2D.markers
Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)
FitFunc=np.polyfit(Times[11:], AvProfVars3[11:, 1], 2, full=False)
Fit = FitFunc[0]*Times[11:]**2 + FitFunc[1]*Times[11:] + FitFunc[2]
dhdt =1.0*(2*FitFunc[0]*Times[11:] + FitFunc[1])/3600
#Fit = FitFunc[0]*Times[120:]**3 + FitFunc[1]*Times[120:]**2 + FitFunc[2]*Times[120:] + FitFunc[3]
#dhdt =1.0*(3*FitFunc[0]*Times[120:]**2 + 2*FitFunc[1]*Times[120:] + FitFunc[2])/3600
deltah = np.subtract(AvProfVars5[:, 2], AvProfVars5[:, 0])
deltah = np.divide(deltah, AvProfVars5[:, 1])
tau = 1.0*rinovals5[:,4]/3600
scaled_time = np.divide(Times, tau)
scaled_dhdt = np.divide(dhdt, rinovals5[11:, 2])
dhdtinvriplt = np.vstack((rinovals5[11:, 1], scaled_dhdt))
dhdtinvriplt = np.transpose(np.vstack((dhdtinvriplt,deltah[11:])))
np.savetxt('/tera/phil/nchaparr/python/Plotting/Mar52014/data/dhdtinvriplt.txt', dhdtinvriplt, delimiter=' ')

#Ax3.plot(dhdtplot[:, 0], dhdtplot[:,1], 'kv', label = '100/10')
#Ax3.plot(dhdtplot0[:, 0], dhdtplot0[:,1], 'ro', label = '100/5')
#Ax3.plot(dhdtplot1[:, 0], dhdtplot1[:,1], 'bo', label = '60/5')
#Ax3.plot(dhdtplot2[:, 0], dhdtplot2[:,1], 'yo', label = '60/2.5')
#Ax3.plot(dhdtplot3[:, 0], dhdtplot3[:,1], 'go', label = '150/5')
#Ax3.plot(dhdtplot4[:, 0], dhdtplot4[:,1], 'co', label = '60/10')
#Ax3.plot(dhdtplot5[:, 0], dhdtplot5[:,1], 'ro', label = '150/10')

#Ax3.plot(Times[11:], HistVars[11:, 1], 'b*', label="h Dist")
#Ax3.plot(rinovals[:,1], deltah[:], 'ko', label = '100/10')
#print len(Times0[7:]), Deltah00[:].shape
Ax3.plot(Times[11:], AvProfVars[11:,1], linestyle='none', marker = 'v', color = 'k', label = '100/10')
Ax3.plot(Times0[7:], AvProfVars0[7:,1], linestyle='none', marker = 'o', color = 'k', label = '100/5')
Ax3.plot(Times[11:], AvProfVars1[11:,1], linestyle='none', marker = 'o', color = 'y', label = '60/5')
Ax3.plot(Times[11:], AvProfVars2[11:,1], linestyle='none', marker = '*', color = 'y', label = '60/2.5')
Ax3.plot(Times[11:], AvProfVars3[11:,1], linestyle='none', marker = 'o', color = 'r', label = '150/5')
Ax3.plot(Times[11:], AvProfVars4[11:,1], linestyle='none', marker = 'v', color = 'y', label = '60/10')
Ax3.plot(Times[11:], AvProfVars5[11:,1], linestyle='none', marker = 'v', color = 'r', label = '150/10')
#Ax3.plot(np.arange(0, .1, .01)[2:10], .20833*np.arange(0, .1, .01)[2:10], 'k--')
#Ax3.plot(np.arange(0, .1, .01)[2:10], np.arange(0, .1, .01)[2:10]**(3.0/2), 'k--')
#Ax3.plot(Times[11:], Fit, 'b-', label="2nd Order Polyfit")
#Ax3.text(6, 1500, r'$ \frac{dh}{dt}  =  %.3f \frac{m}{s} $' %(1.0*M/3600),  fontdict=None, withdash=False, fontsize = 15)
#Ax3.set_ylim(0, 2500)
Ax3.legend(loc = 'upper left', prop={'size': 10})
#Ax3.set_title(r'$\Delta h (Flux)\ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$Scaled \ Time \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\frac{\Delta h}{h} \ vs \ Ri^{-1}$', fontsize=20)
#Ax3.set_title(r'$Ri^{-1} \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\Delta \theta \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\overline{\theta} \ vs \ Time$', fontsize=20)
#Ax3.set_xlabel(r"$\frac{Time}{\tau}$", fontsize=20)
#Ax3.set_ylabel(r"$\frac{\Delta h}{h}$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta h (m)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta \theta (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\overline{ \theta} (K)$", fontsize=20)
Ax3.set_ylabel(r"$h \ (m)$", fontsize=20)
Ax3.set_xlabel(r"$Time (hrs)$", fontsize=20)
#Ax3.set_ylabel(r"$h \ (m)$", fontsize=20)
#plt.ylim(0, .61)
plt.show()





    
    
