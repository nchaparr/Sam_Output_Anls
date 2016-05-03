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
        Deltah = points.Deltah_over_h()
        HistVars = points.HistVars()
        AvProfVars = points.AvProfVars() #chage to AvProfVars_old() for old limits
        gm_vars = points.gm_vars()

        if Run_Date_List[i]=="Nov302013":
            print "Nov302013", gm_vars[11:,3].shape, AvProfVars[11:,8].shape, AvProfVars[11:, 3].shape
            Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 3], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 4], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 1], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 0], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times0[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 4], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times0[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 5], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times0[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 6], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times0[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 7], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
        else:
            print gm_vars[11:,3].shape, AvProfVars[11:,8].shape, AvProfVars[11:, 3].shape
            Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 3], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 4], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 1], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(AvProfVars[11:,8], np.divide(AvProfVars[11:, 0], gm_vars[11:, 3]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 4], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 5], gm_vars[11:, 0]), legend_list[i], label = label_list[i])
            #Ax3.plot(np.multiply(Times[11:],gm_vars[11:,1]), np.divide(AvProfVars[11:, 7], gm_vars[11:, 0]), legend_list[i], label = label_list[i])

#Ax3.text(33.5, 0.8, r"$z_{g0}$", size=20)
#Ax3.text(33.5, .95, r"$z_{f0}$", size=20)
#Ax3.text(33.5, 1.1, r"$z_{f}$", size=20)
#Ax3.text(33.5, 1.15, r"$z_{g}$", size=20)
                    
#Ax3.legend(loc = 'lower right', prop={'size': 10}, numpoints=1)
#Ax3.set_title(r'GM Height vs Time Scaling', fontsize=20)
Ax3.set_ylabel(r"$c_{\delta}$", fontsize=20)
Ax3.set_xlabel(r"$z_{enc} / L_{0}}$", fontsize=20)
plt.ylim(0.6, 1.4)
plt.show()





    
    
