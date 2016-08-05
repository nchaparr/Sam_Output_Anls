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

#Main Part -- pulling points and plotting them
label_list = ['100/10', '100/5', '60/5', '60/2.5', '150/5', '60/10', '150/10']
legend_list = ['kv', 'ko', 'yo', 'y*', 'ro', 'yv', 'rv']
Run_Date_List = ["Dec142013", "Nov302013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]


for i in range(len(label_list)):
    if i<99:
        points = For_Plots(Run_Date_List[i])
        rinovals = points.rinovals()
        Deltah = points.Deltah(8, 1)#upper, lower, scale
        HistVars = points.HistVars()
        AvProfVars = points.AvProfVars() #chage to AvProfVars_old() for old limitsx
        gm_vars = points.gm_vars()


        #Jan152014_1: 16:21=np.nan, plot11:29
        #Mar12014; 12:17=nan
        if Run_Date_List[i]=="Jan152014_1":
            #gm_vars[16:21, 3]=np.nan
            #gm_vars[16:21:, 0]=np.nan
            #AvProfVars[16:21, 8]=np.nan
            #AvProfVars[16:21, 1]=np.nan
            #AvProfVars[16:21, 2]=np.nan
            AvProfVars[14:17, 0]=np.nan
            #print "Nov302013", gm_vars[11:,3].shape, AvProfVars[11:,8].shape, AvProfVars[11:, 3].shape
            #Ax3.plot(np.divide(gm_vars[11:29, 3], gm_vars[11:29, 0]), rinovals[11:29, 1], legend_list[i], label = label_list[i])
            
            Ax3.plot(np.divide(gm_vars[11:29, 3], gm_vars[11:29, 0]), np.divide(AvProfVars[11:29, 2], gm_vars[11:29, 3]), legend_list[i], label = label_list[i])
            Ax3.plot(np.divide(gm_vars[11:29, 3], gm_vars[11:29, 0]), np.divide(AvProfVars[11:29, 4], gm_vars[11:29, 3]), legend_list[i], label = label_list[i])
            #Ax3.loglog(np.multiply(rinovals[11:29, 10]**(0.5),np.divide(gm_vars[11:29, 3], gm_vars[11:29, 0])), np.divide(Deltah[11:29], gm_vars[11:29,3]), legend_list[i], label = label_list[i], markersize=10)
            #Ax3.plot(np.divide(gm_vars[11:29,3], gm_vars[11:29,0]),rinovals[11:29, 10],  legend_list[i], label = label_list[i], markersize=10)
        elif Run_Date_List[i]=="Mar12014":
            gm_vars[12:17, 3]=np.nan
            gm_vars[12:17:, 0]=np.nan
            AvProfVars[12:17, 8]=np.nan
            AvProfVars[12:17, 1]=np.nan
            AvProfVars[12:17, 2]=np.nan
            AvProfVars[12:17, 0]=np.nan
            #print "Nov302013", gm_vars[11:,3].shape, AvProfVars[11:,8].shape, AvProfVars[11:, 3].shape
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 3], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), rinovals[11:, 1], legend_list[i], label = label_list[i])
            Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 2], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 3], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 5], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.loglog(np.multiply(rinovals[11:29, 10]**(0.5), np.divide(gm_vars[11:29,3], gm_vars[11:29,0])), np.divide(Deltah[11:29],gm_vars[11:29,3]), legend_list[i], label = label_list[i], markersize=10)
            #Ax3.plot(np.divide(gm_vars[11:29,3], gm_vars[11:29,0]),rinovals[11:29, 10],  legend_list[i], label = label_list[i], markersize=10)

        elif Run_Date_List[i]=="Mar52014":
            #gm_vars[12:17, 3]=np.nan
            #gm_vars[12:17:, 0]=np.nan
            #AvProfVars[12:17, 8]=np.nan
            #AvProfVars[12:17, 1]=np.nan
            #AvProfVars[12:17, 2]=np.nan
            AvProfVars[28:32, 0]=np.nan
            #print "Nov302013", gm_vars[11:,3].shape, AvProfVars[11:,8].shape, AvProfVars[11:, 3].shape
            Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 3], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 2], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), rinovals[11:, 1], legend_list[i], label = label_list[i])
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 3], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 5], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.loglog(np.multiply(rinovals[11:, 10]**(0.5),np.divide(gm_vars[11:, 3], gm_vars[11:, 0])), np.divide(Deltah[11:],gm_vars[11:,3]), legend_list[i], label = label_list[i], markersize=10)
            #Ax3.plot(np.divide(gm_vars[11:,3], gm_vars[11:,0]),rinovals[11:, 10],  legend_list[i], label = label_list[i], markersize=10)         
        elif Run_Date_List[i]=="Nov302013":
            AvProfVars[13, 8]=np.nan
            AvProfVars[13, 1]=np.nan
            AvProfVars[13, 2]=np.nan
            AvProfVars[13, 0]=np.nan
            AvProfVars[15:17, 8]=np.nan
            AvProfVars[15:17, 1]=np.nan
            AvProfVars[15:17, 2]=np.nan
            AvProfVars[15:17, 0]=np.nan
            AvProfVars[24:26, 8]=np.nan
            AvProfVars[24:26, 1]=np.nan
            AvProfVars[24:26, 2]=np.nan
            AvProfVars[24:26, 0]=np.nan            
            #print "Nov302013", gm_vars[11:,3].shape, AvProfVars[11:,8].shape, AvProfVars[11:, 3].shape
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 3], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), rinovals[11:, 1], legend_list[i], label = label_list[i])
            Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 2], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 3], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 5], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 4], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 1], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 0], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times0[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 4], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times0[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 5], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times0[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 6], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times0[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 7], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.loglog(np.multiply(rinovals[7:, 10]**(0.5),np.divide(gm_vars[7:, 3], gm_vars[7:, 0])), np.divide(Deltah[7:],gm_vars[7:,3]), legend_list[i], label = label_list[i], markersize=10)
            #Ax3.plot(np.divide(gm_vars[7:,3], gm_vars[7:,0]),rinovals[7:, 10],  legend_list[i], label = label_list[i], markersize=10)         
        else:
            print(Run_Date_List[i], label_list[i])
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 3], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), rinovals[11:, 1], legend_list[i], label = label_list[i])
            Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 2], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 4], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.divide(gm_vars[11:, 3], gm_vars[11:, 0]), np.divide(AvProfVars[11:, 5], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 4], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 1], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 0], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 4], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 5], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 7], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.loglog(np.multiply(rinovals[11:, 10]**(0.5), np.divide(gm_vars[11:, 3], gm_vars[11:, 0])), np.divide(Deltah[11:],gm_vars[11:,3]), legend_list[i], label = label_list[i], markersize=10)
            #Ax3.plot(np.divide(gm_vars[11:,3], gm_vars[11:,0]),rinovals[11:, 10],  legend_list[i], label = label_list[i], markersize=10)


#xes = np.arange(6, 25, 1)
#x1es = np.arange(.029, .04, .0001)
#ys = 2.5*xes**(.5)
#ys= (.4**.5)*xes**(-2.0/3)
#Ax3.loglog(xes, ys, 'k--')
#Ax3.loglog(x1es, ys1, 'k--')


#Ax3.text(32.5, 1.25, r"$z_{g1}$", size=20)
Ax3.text(32.5, 1.2, r"$z_{f1}$", size=20)
Ax3.text(32.5, 1.2, r"$z_{g1GM}$", size=20)
Ax3.text(32.5, .95, r"$z_{f0}$", size=20)
#Ax3.text(32.5, 1.1, r"$z_{f}$", size=20)
#Ax3.text(32.5, 1.15, r"$z_{g}$", size=20)
#Ax3.text(32.5, 0.87, r"$z_{g0}$", size=20)

Ax3.text(6, 1.4, '(a)', fontsize=30)


#Ax3.text(6.5, .2, r'$y = \sqrt{0.4}x^{-\frac{2}{3}}$',  fontdict=None, withdash=False, fontsize = 25, rotation=-8)
   
#Ax3.set_xlabel(r"$c_{\delta}(z_{enc}/L_{0})$", fontsize=30)                   
#Ax3.legend(loc = 'lower right', prop={'size': 10}, numpoints=1)
#Ax3.set_title(r'GM Height vs Time Scaling', fontsize=20)
#Ax3.set_ylabel(r"$\delta/z_{enc}$", fontsize=30)
Ax3.set_ylabel(r"$z / z_{enc}$", fontsize=30)
Ax3.set_xlabel(r"$z_{enc} / L_{0}}$", fontsize=30)
Ax3.tick_params(axis="both", labelsize=20)
plt.ylim(0.6, 1.5)
plt.xlim(5, 35)
plt.tight_layout()
plt.show()

    
    
