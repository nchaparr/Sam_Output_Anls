import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
from nchap_class import *
import nchap_fun as nc
from matplotlib.lines import Line2D
from scipy.stats import linregress
from matplotlib import rcParams
rcParams.update({'font.size': 10})

"""
   plots entrainment rate plot
  
"""

dump_time_list, Times = Make_Timelists(1, 600, 28800)
Times = np.array(Times)
dump_time_list0, Times0 = Make_Timelists(1, 900, 28800)
Times0 = np.array(Times0)
Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)


#Main Part -- pulling points and plotting them
label_list = ['100/10', '100/5', '60/5', '60/2.5', '150/5', '60/10', '150/10']
legend_list = ['kv', 'ko', 'yo', 'y*', 'ro', 'yv', 'rv']
Run_Date_List = ["Dec142013", "Nov302013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]

regress_x=np.array([0])
regress_y=np.array([0])
for i in range(len(label_list)):
    if i<99:
        points = For_Plots(Run_Date_List[i])
        rinovals = points.rinovals()
        Deltah = points.Deltahf_over_zf()
        HistVars = points.HistVars()
        AvProfVars = points.AvProfVars()
        
        if Run_Date_List[i] == "Nov302013":
             points.Get_and_save_dhdt(Times0[7:], AvProfVars[7:, 1], rinovals[7:, 2], rinovals[7:, 1])
             scaled_we_plot = points.scaled_we_plot()
             Ax3.loglog(scaled_we_plot[0, :], scaled_we_plot[1, :], legend_list[i], label = label_list[i], markersize=12)
             regress_x=np.hstack((regress_x,scaled_we_plot[0, :]))
             regress_y=np.hstack((regress_y,scaled_we_plot[1, :]))             
        
        elif Run_Date_List[i] == "Jan152014_1":
             points.Get_and_save_dhdt(Times[11:29], AvProfVars[11:29, 1], rinovals[11:29, 2], rinovals[11:29, 1])
             scaled_we_plot = points.scaled_we_plot()
             Ax3.loglog(scaled_we_plot[0, :], scaled_we_plot[1, :], legend_list[i], label = label_list[i], markersize=12)
             regress_x=np.hstack((regress_x,scaled_we_plot[0, :]))
             regress_y=np.hstack((regress_y, scaled_we_plot[1, :]))
        
        elif Run_Date_List[i] == "Mar12014":
             points.Get_and_save_dhdt(Times[11:], AvProfVars[11:, 1], rinovals[11:, 2], rinovals[11:, 1])
             scaled_we_plot = points.scaled_we_plot()
             Ax3.loglog(scaled_we_plot[0, :], scaled_we_plot[1, :], legend_list[i], label = label_list[i], markersize=12)            
             regress_x=np.hstack((regress_x, scaled_we_plot[0, :]))
             regress_y=np.hstack((regress_y, scaled_we_plot[1, :]))
        
        else:
             points.Get_and_save_dhdt(Times[11:], AvProfVars[11:, 1], rinovals[11:, 2], rinovals[11:, 1])
             scaled_we_plot = points.scaled_we_plot()
             Ax3.loglog(scaled_we_plot[0, :], scaled_we_plot[1, :], legend_list[i], label = label_list[i], markersize=12)
             regress_x=np.hstack((regress_x, scaled_we_plot[0, :]))
             regress_y=np.hstack((regress_y, scaled_we_plot[1, :]))

regress_x=np.delete(regress_x,0)
regress_y=np.delete(regress_y,0)

#log10_regress_x=np.log10(regress_x)
#log10_regress_y=np.log10(regress_y)
#log10_regress_results = linregress(log10_regress_x, log10_regress_y)

regress_results = linregress(regress_x, regress_y)

xs = np.arange(0.02, .09, .0001)
#ys=  log10_regress_results.intercept + log10_regress_results.slope*xs
ys = regress_results.intercept + regress_results.slope*xs

Ax3.loglog(xs, ys, 'k--')

Ax3.text(.03, .005, '(b)',  fontdict=None, withdash=False, fontsize = 30)
Ax3.set_ylabel(r"$\frac{w_{e}}{w^{*}}$", fontsize=40)
Ax3.set_xlabel(r"$Ri_{\Delta g}^{-1}$", fontsize=35)
plt.xlim(.028, .12)
plt.ylim(0.0045, 0.05)
Ax3.set_yticks([.005, .01, .05])
Ax3.set_yticklabels([.005, .01, .05])
Ax3.set_xticks([.03, .06, .12])
Ax3.set_xticklabels([.03, .06, .12])
Ax3.tick_params(axis="both", labelsize=25)
plt.tight_layout()
plt.show()





    
    
