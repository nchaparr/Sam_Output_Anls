import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import Make_Timelists
from nchap_class import For_Plots
import nchap_fun as nc
from matplotlib.lines import Line2D

from matplotlib import rcParams
rcParams.update({'font.size': 10})

"""
   plots GM eq 26
  
"""

dump_time_list, Times = Make_Timelists(1, 600, 28800)
Times = np.array(Times)
dump_time_list0, Times0 = Make_Timelists(1, 900, 28800)
Times0 = np.array(Times0)

#plot

Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)


#Main Part -- pulling points and plotting them
label_list = ['100/10', '100/5', '60/5', '60/2.5', '150/5', '60/10', '150/10']
legend_list = ['kv', 'ko', 'yo', 'y*', 'ro', 'yv', 'rv']
Run_Date_List = ["Dec142013", "Nov302013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]


for i in range(len(label_list)):
    if i<99:
        
        #pull variables from text files 
        points = For_Plots(Run_Date_List[i])
        rinovals = points.rinovals()
        gm_vars = points.gm_vars()
        heights = points.AvProfVars()
        
        #name and calculate points for plot
        z1_GM = heights[:,8]
        z_g=heights[:,1] # [elbot_dthetadz, h, eltop_dthetadz, elbot_flux, h_flux, eltop_flux, deltatheta, mltheta, z1_GM]
        delta = z1_GM - z_g
        z_enc=gm_vars[:, 3]
        L0 = gm_vars[:, 0]
        c_squared_delta = rinovals[:, 10]
        c_delta = c_squared_delta**.5

        Xs =  c_delta*(z_enc/L0) 
        Ys = delta/z_enc
        #Ys = points.Deltah_over_h(8, 1, 7)

        if Run_Date_List[i] == "Nov302013":
             #clean up
             Ys[13] = np.nan
             Ys[15:17] = np.nan
             Ys[24:26] = np.nan
             #plot 
             Ax3.loglog(Xs[7:], Ys[7:], legend_list[i], label = label_list[i], markersize=10)
             
        elif Run_Date_List[i] == "Jan152014_1":
             #clean up
             Ys[16:21] = np.nan
             #plot
             Ax3.loglog(Xs[11:29], Ys[11:29], legend_list[i], label = label_list[i], markersize=10)
        
        elif Run_Date_List[i] == "Mar12014":
             #clean up
             Ys[11:17] = np.nan
             #plot
             Ax3.loglog(Xs[11:29], Ys[11:29], legend_list[i], label = label_list[i], markersize=10)

        else:     
             Ax3.loglog(Xs[11:], Ys[11:], legend_list[i], label = label_list[i], markersize=10)

xs = np.arange(6, 25, 1)
ys= (.4**.5)*xs**(-2.0/3)
Ax3.loglog(xs, ys, 'k--')

Ax3.text(6.6, .2, r'$y = \sqrt{0.4}x^{-\frac{2}{3}}$',  fontdict=None, withdash=False, fontsize = 25, rotation=-8)
Ax3.set_xlabel(r"$c_{\delta}(z_{enc}/L_{0})$", fontsize=20)
Ax3.set_ylabel(r"$\delta/z_{enc}$", fontsize=20)
#plt.ylim(0, 0.2)
plt.xlim(0, 22)
Ax3.tick_params(axis="both", labelsize=20)
plt.tight_layout()
plt.show()





    
    
