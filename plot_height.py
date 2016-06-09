import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import Make_Timelists
from nchap_class import For_Plots
#from Sam_Output_Anls.Make_Timelist import Make_Timelists
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
#from Sam_Output_Anls.nchap_class import For_Plots
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
        rinovals = points.rinovals() #rinovals_old
        gm_vars = points.gm_vars()
        Deltah = points.Deltah_over_h(2, 0, 1) #change for deltahinvri plots [elbot_dthetadz, h, eltop_dthetadz, elbot_flux, h_flux, eltop_flux, deltatheta, mltheta, z1_GM] Deltah_old
        HistVars = points.HistVars()  
        AvProfVars = points.AvProfVars()  # AvProfVars_old()
        if Run_Date_List[i] == "Nov302013":
             Deltah[13] = np.nan
             Deltah[15:17] = np.nan
             Deltah[24:26] = np.nan
             #Ax3.plot(np.divide(gm_vars[7:, 3], gm_vars[7:, 0]), rinovals[7:,10], legend_list[i], label = label_list[i], markersize=10) 
             #Ax3.loglog(np.multiply(rinovals[7:, 10]**(0.5),np.divide(gm_vars[7:, 3], gm_vars[7:, 0])), Deltah[7:], legend_list[i], label = label_list[i], markersize=10) #for GM comparison plot
             Ax3.loglog(rinovals[7:, 1], Deltah[7:], legend_list[i], label = label_list[i], markersize=10) #for our invri deltah plots
             #Ax3.loglog(np.divide(gm_vars[7:, 3], gm_vars[7:, 0]), Deltah[7:], legend_list[i], label = label_list[i], markersize=10)
        elif Run_Date_List[i] == "Jan152014_1":
    #TODO: alternative starting index for Nov302013
             Deltah[16:21] = np.nan
             #print Deltah
             #Ax3.plot(np.divide(gm_vars[11:29,3], gm_vars[11:29,0]), rinovals[11:29,10], legend_list[i], label = label_list[i], markersize=10)
             #Ax3.loglog(np.multiply(rinovals[11:29, 10]**(0.5), np.divide(gm_vars[11:29,3], gm_vars[11:29,0])), Deltah[11:29], legend_list[i], label = label_list[i], markersize=10) #for GM comparison plot
             Ax3.loglog(rinovals[11:29, 1], Deltah[11:29], legend_list[i], label = label_list[i], markersize=10) #for our deltah invrino plot
             #Ax3.loglog(np.divide(gm_vars[11:29,3], gm_vars[11:29,0]), Deltah[11:29], legend_list[i], label = label_list[i], markersize=10)
        elif Run_Date_List[i] == "Mar12014":
    #TODO: alternative starting index for Nov302013
             Deltah[11:17] = np.nan
             #print Deltah
             #Ax3.loglog(np.multiply(rinovals[11:29, 10]**(0.5), np.divide(gm_vars[11:29, 3], gm_vars[11:29, 0])), Deltah[11:29], legend_list[i], label = label_list[i], markersize=10) #for GM comparison plot
             Ax3.loglog(rinovals[11:29, 1], Deltah[11:29], legend_list[i], label = label_list[i], markersize=10) #for our deltah invrino plot
             #Ax3.loglog(np.divide(gm_vars[11:29, 3], gm_vars[11:29, 0]), Deltah[11:29], legend_list[i], label = label_list[i], markersize=10)
             #Ax3.plot(np.divide(gm_vars[11:29, 3], gm_vars[11:29, 0]), rinovals[11:29,10], legend_list[i], label = label_list[i], markersize=10)
        else:     
             #print  rinovals[11:, 10]
             zenc_over_L0 = np.divide(gm_vars[11:, 3], gm_vars[11:, 0])
             #Ax3.loglog(np.multiply(rinovals[11:, 10]**(0.5), zenc_over_L0), Deltah[11:], legend_list[i], label = label_list[i], markersize=10) #for GM comparison plot
             Ax3.loglog(rinovals[11:, 1], Deltah[11:], legend_list[i], label = label_list[i], markersize=10) #for our invro deltah plot
             #Ax3.loglog(zenc_over_L0, Deltah[11:], legend_list[i], label = label_list[i], markersize=10)
             #Ax3.plot(zenc_over_L0, rinovals[11:,10], legend_list[i], label = label_list[i], markersize=10)

#xes = np.arange(6, 25, 1) #for gm comparison plot
#ys= (.4**.5)*xes**(-2.0/3)#for gm comparison plot
#Ax3.loglog(xes, ys, 'k--')#for gm comparison plot

xes = np.arange(.04, .075, .0001)#for our deltah invri plot
x1es = np.arange(.029, .04, .0001)#for our deltah invri plot
ys = 2.5*xes**(.5)#for our deltah invri plot
ys1= 12.5*x1es**(1)#for our deltah invri plot
Ax3.loglog(xes, ys, 'k--')#for our deltah invri plot
Ax3.loglog(x1es, ys1, 'k--')#for our deltah invri plot

#Ax3.plot(np.arange(0, .1, .01)[2:10], np.arange(0, .1, .01)[2:10]**(3.0/2), 'k--')
#Ax3.plot(Times[11:], Fit, 'b-', label="2nd Order Polyfit")
#Ax3.text(6.6, .2, r'$y = \sqrt{0.4}x^{-\frac{2}{3}}$',  fontdict=None, withdash=False, fontsize = 25, rotation=-8) #for GM comparison plot
Ax3.text(.028, .49, r'$b = -1$',  fontdict=None, withdash=False, fontsize = 25, rotation=35)#for our deltah invri plots
Ax3.text(.048, .67, r'$b = -\frac{1}{2}$',  fontdict=None, withdash=False, fontsize = 25, rotation=22)

#Ax3.text(.085, .88, r'(d)', fontsize=30)
#Ax3.set_ylim(0, 2500)
#Ax3.legend(loc = 'lower right', prop={'size': 14}, numpoints=1)
#Ax3.set_title(r'$\Delta h (Flux)\ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$Scaled \ Time \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\frac{\Delta h}{h} \ vs \ Ri^{-1}$', fontsize=20)
#Ax3.set_title(r'$Ri^{-1} \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\Delta \theta \ vs \ Time$', fontsize=20)
#Ax3.set_title(r'$\overline{\theta} \ vs \ Time$', fontsize=20)
#Ax3.set_xlabel(r"$\frac{Time}{\tau}$", fontsize=20)
#Ax3.set_ylabel(r"$\frac{z}{h}$", fontsize=20)
#Ax3.set_ylabel(r"$\frac{\Delta z}{z_g}$", fontsize=30)
#Ax3.set_ylabel(r"$\frac{w_{e}}{w^{*}}$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta h (m)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta \theta (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\overline{ \theta} (K)$", fontsize=20)
#Ax3.set_ylabel(r"$\Delta h \ (m)$", fontsize=20)
#Ax3.set_xlabel(r"$c_{\delta}(z_{enc}/L_{0})$", fontsize=20)#gm comp plot
#Ax3.set_xlabel(r"$Ri_{\Delta}^{-1}$", fontsize=30)
#Ax3.set_xlabel(r"$\gamma \frac{\Delta h}{\Delta \theta}$", fontsize=20)
#Ax3.set_ylabel(r"$\delta/z_{enc}$", fontsize=20) #gm comp plot
#plt.ylim(0, 0.2)
#Ax3.set_yticks([0, 0.1, 0.2])
#Ax3.set_yticklabels([0.2, 0.4, 0.6, 0.8, 1])
#plt.xlim(6, 22) #gm comp plot
#Ax3.set_xticks([0, 10])
#Ax3.set_xticklabels([0.02, 0.04, 0.06, 0.08, .1])

#for deltah invri plots
plt.ylim(0.2, 1)
Ax3.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
Ax3.set_yticklabels([0.2, 0.4, 0.6, 0.8, 1])
plt.xlim(.02, .1)
Ax3.set_xticks([0.02, 0.04, 0.06, 0.08, .1])
Ax3.set_xticklabels([0.02, 0.04, 0.06, 0.08, .1])


Ax3.tick_params(axis="both", labelsize=20)
plt.tight_layout()
plt.show()





    
    
