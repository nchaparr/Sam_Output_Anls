import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import site
site.addsitedir('/tera/phil/nchaparr/python')
import nchap_fun as nc

from Make_Timelist import *

"""
   for plotting EL limits per average profile on Scaling Diagram

"""

#TODO: edit so that it only plots the scaing diagram from av prof values


#create lists of txt file to loop over
dump_time_list, Times = Make_Timelists(1, 600, 28800)

ELLimits = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec252013/data/ELLims_hist")
#ELLimits1 = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Sep302013/data/ELLims_hist1")
  
AvProfVars = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec252013/data/AvProfLims")
rinovals = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec252013/data/invrinos")

#stull_data = np.genfromtxt('/tera/phil/nchaparr/python/Plotting/July1112013/data/stull_vars.txt')
#mol_etc = nc.from_lmo()
#hoverL = -np.divide(AvProfVars[:, 1], mol_etc[:, 0])
#print mol_etc
#plot the heights vs time
Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)
#Ax3.plot(rinovals[60:, 1], np.divide(AvProfVars[:,0], AvProfVars[:,1])[60:],'yo', label = 'nchap')
#Ax3.plot(rinovals[60:, 1], np.divide(AvProfVars[:,2], AvProfVars[:,1])[60:], 'go')
#Ax3.plot(-stull_data[2, :], stull_data[0, :], 'y*', label = 'from stull')
#Ax3.plot(-stull_data[2, :], stull_data[1, :], 'g*')
#Ax3.plot(-stull_data[2, 0:6], np.zeros_like(stull_data[2, 0:6])+1.2, 'k--', label = 'Holt Diagram')
#Ax3.plot(-stull_data[2, 0:6], np.zeros_like(stull_data[2, 0:6])+.8, 'k--')
Ax3.plot(Times[11:], np.divide(AvProfVars[:,0], AvProfVars[:,1])[11:],'ko', label = r"$from \ \frac{\partial \overline{\theta}}{\partial z}$")
Ax3.plot(Times[11:], np.divide(AvProfVars[:,2], AvProfVars[:,1])[11:], 'ko')
#Ax3.plot(Times[11:], np.divide(AvProfVars[:,3], AvProfVars[:,1])[11:],'b*', label = r"$from \ \overline{w^{'}\theta^{'}}$")
#Ax3.plot(Times[11:], np.divide(AvProfVars[:,5], AvProfVars[:,1])[11:], 'b*')
Ax3.plot(Times[11:], np.divide(ELLimits[:,0], AvProfVars[:,1])[11:],'b*', label = r"$from \ Percentiles$")
Ax3.plot(Times[11:], np.divide(ELLimits[:,2], AvProfVars[:,1])[11:], 'b*')

plt.legend(loc = 'lower right', prop={'size':8})
Ax3.set_title(r"$Scaled \ EL \ Limits$", fontsize=20)
Ax3.set_xlabel(r"$Time \ (hrs)$", fontsize=20)
Ax3.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
plt.ylim(0, 1.5)
plt.show()



    
    
