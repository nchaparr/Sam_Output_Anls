from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
import matplotlib.pyplot as plt
import site
from datetime import datetime
#site.addsitedir('/tera/phil/nchaparr/python')
#from nchap_fun import *
from Make_Timelist import *
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from nchap_class import *

"""
    plotting output from an ensemble of runs
    writes out text files for get_dtheta and get_rhino
      
"""

def Main_Fun(dump_time):
     
     """Pulls output from an ensemble cases, gets ensemble averages and perturbations and
     their horizontal averages

    Arguments:
    dump_time -- time of output eg '0000000720'

    Returns:
    var_bar -- array of horizontally averaged, ensemble averages or perturbations (covariances)
    height -- array of height values
    
    """
     date = "Dec142013" #TODO: this should be an argument passed to Main_Fun
     #pulls data using class Get_Var_Arrays1     
     Vars = Get_Var_Arrays1("/tera2/nchaparr/"+date+"/runs/sam_case", "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_", dump_time)
     thetas_list, press_list = Vars.get_thetas() 
     wvels_list= Vars.get_uvelperts()
     height = Vars.get_height()
                  
     #get arrays of ensemble averaged variables
     ens_avthetas = nc.Ensemble1_Average(thetas_list) #TODO: 1 could be removed from this in nchap_fun
     #avthetas = nc.Horizontal_Average(thetas_list[0])  #for fake perturbations from a single case
     #ens_avthetas = np.zeros_like(theta)
     #for i in range(312):
     #     ens_avthetas[i,:,:] = avthetas[i]
     #ens_avwvels = nc.Ensemble1_Average(wvels_list)
     ens_press = nc.Ensemble1_Average(press_list)
              
     #now get the perturbations
     #wvelperts_list = []
     wvelpertsq_list = Vars.get_sqvel('v')
     #thetaperts_list = []
     wvelthetaperts_list = Vars.get_wvelthetaperts()           
               
     #and ensemble average them
     #print (datetime.now()-startTime)
     
     #ens_avthetaperts = nc.Ensemble1_Average(thetaperts_list)
     #ens_avwvelperts = nc.Ensemble1_Average(wvelperts_list)
     ens_avwvelpertsq = nc.Ensemble1_Average(wvelpertsq_list)     
     rtens_avwvelpertsq = np.sqrt(ens_avwvelpertsq)
     ens_avwvelthetaperts = nc.Ensemble1_Average(wvelthetaperts_list)
     
     #loop over ensemple averages and horizontally average
     #wvel_bar = nc.Horizontal_Average(ens_avwvels)     
     theta_bar = nc.Horizontal_Average(ens_avthetas)
     #wvelpert_bar = Horizontal_Average(ens_avwvelperts)
     wvelthetapert_bar = nc.Horizontal_Average(ens_avwvelthetaperts)
     
     #for testing w scale -- slows script down
     rtwvelpertsq_bar = nc.Horizontal_Average(rtens_avwvelpertsq)
     
     #save text files, TODO: make more purdy
     #np.savetxt('/tera/phil/nchaparr/python/Plotting/"+date+"/data/flux_quads' + dump_time, np.transpose(np.array([upwarm_bar, downwarm_bar, upcold_bar, downcold_bar])), delimiter=' ')
     #np.savetxt('/tera/phil/nchaparr/python/Plotting/'+date+'/data/wvelthetapert'+dump_time, wvelthetapert_bar, delimiter=' ')
     #np.savetxt('/tera/phil/nchaparr/python/Plotting/'+date+'/data/theta_bar'+dump_time, theta_bar, delimiter=' ')
     #np.savetxt('/tera/phil/nchaparr/python/Plotting/'+date+'/data/heights'+dump_time, height, delimiter=' ')
     #np.savetxt('/tera/phil/nchaparr/python/Plotting/'+date+'/data/press'+dump_time, ens_press, delimiter=' ')
     #np.savetxt('/tera/phil/nchaparr/python/Plotting/"+date+"/data/tracers'+dump_time, tracer_bar, delimiter=' ')
     #np.savetxt('/tera/phil/nchaparr/python/Plotting/'+date+'/data/rootmeanvsq'+dump_time, rtwvelpertsq_bar, delimiter=' ')
     
     return rtwvelpertsq_bar, height


go_ahead = np.int(raw_input('have you changed the write out folder paths? 1 or 0: '))
if go_ahead == 1:

     #MLZero_Vars = np.genfromtxt('/tera/phil/nchaparr/python/Plotting/Nov42013/data/MLZero_Vars.txt')
     
     dump_time_list, Times = Make_Timelists(1, 600, 28800)
     

     #set up plot
     theAx = nc.Do_Plot(3, r"Horizontally Averaged, Ensemble Averaged $w^{'}\theta^{'}$ Profiles", 'Height (m)', r"$\overline{w^{'}\theta^{'}}$ (mK/s)", 111)
     #get horizontally averaged ensemble averaged variable and plot
     colorlist=['k', 'b', 'c', 'g', 'r', 'm', 'y', '.75']
     for i in range(48):
          if np.mod(i+1, 1)==0:
               print i
               
               #make plots for MLZero
               color = colorlist[int(1.0*i/6)]
               
               #theAx.plot([MLZero_Vars[i, 0], MLZero_Vars[i, 0]], [0, MLZero_Vars[i, 1]], color+'--')
               #theAx.plot([MLZero_Vars[i, 0], MLZero_Vars[i, 0]+MLZero_Vars[i, 2]], [MLZero_Vars[i, 1], MLZero_Vars[i, 1]], color+'--')
               #theAx.plot([MLZero_Vars[0:i-2]+MLZero_Vars[2:i-2], ], [MLZero_Vars[1:i-2], MLZero_Vars[1:i-2]], '--')
               wvelpert_bar, height = Main_Fun(dump_time_list[i])
               print wvelpert_bar.shape, height.shape
               #wvelpert_bar[0] = np.nan #for getting rid of boundary value effect
               #wvelpert_bar = np.multiply(wvelpert_bar, np.zeros_like(wvelpert_bar)+1004) #for converting to watts/m2
               
               theAx.plot(wvelpert_bar, height, color, label = str(Times[i]) + 'hrs')
               #Encroachment/Thermodynamic Model #TODO: delete this? when sure no longer needed
               #X = 1.0*(i+120)*60*80/(1004*1.225)#Wm-2s
               #h = np.sqrt(200*X)
               #T = 1.0*(h + 60000)/200
               #print 'Encroachment Values', h, T, 1.0*(i+120)*60*60/(1004*1.225), [T, 0], [T, h]
               #theAx.plot([T, T], [0, h], ':')
               #theAx.text(300, 1500, '-- Zero Order Model',  fontdict=None, withdash=False)
               #theAx.text(300, 1400, ': Encroachment',  fontdict=None, withdash=False)
     zeros = np.zeros(len(height))
     theAx.plot(zeros, height)

     #get initial sounding data ie potential temperature
     array = np.genfromtxt('/tera/phil/nchaparr/python/Pert_Files/snd1')
     height_0 = array[:, 0]
     theta_0 = array[:, 1]
     #f=interp1d(height_0, theta_0) #not sure i need this

     #Now plot inital sounding
     #top_index = np.where(height <= 1130)[0][-1]
     #theta_0 = f(height[0:top_index])
     #theAx.plot(theta_0, height_0, label = 'initial sounding')
     plt.legend(loc = 'lower right', prop={'size':8})
     plt.ylim(50, 2000)
     plt.xlim(300, 320)
     plt.show()

else:
    print 'need to update write out folders' 
     


    
    
