from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import site
from datetime import datetime
site.addsitedir('/tera/phil/nchaparr/python')
import nchap_fun as nc
from matplotlib.colors import Normalize
from Make_Timelist import *
from nchap_class import *

#site.addsitedir('/tera2/nchaparr/Dcc252013/hist2d')
#from nchap_2dhist import *


"""  
     For plotting Flux quadrants    
"""

def Main_Fun(date, dump_time, hflux):
     
     """Pulls output from an ensemble cases, gets ensemble averages and perturbations and
     their horizontal averages

    Arguments:
    dump_time -- time of output eg '0000000720'

    Returns:
    var_bar -- 64 array of horizontally averaged, ensemble averages or perturbations (covariances)
    
    """
     #create list of filenames for given dump_time     
     ncfile_list = ["/newtera/tera/phil/nchaparr/tera2_cp/nchaparr/"+date+"/runs/sam_case" + str(i+1) + "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_" + dump_time + ".nc" for i in range(10)]

     #create lists for variable arrays from each case
     upwarm_list = []
     downwarm_list = []
     upcold_list = []
     downcold_list = []
     wvelperts_list = []
     thetaperts_list = []
     wvelthetaperts_list = []

     #get velocity perts and thetas
     Vars =  Get_Var_Arrays1("/newtera/tera/phil/nchaparr/tera2_cp/nchaparr/"+date+"/runs/sam_case", "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_", dump_time)         
     thetas_list, press_list = Vars.get_thetas()     
     wvels_list = Vars.get_wvelperts()          
     height = Vars.get_height()
     
     #get arrays of enseble averaged variables
     ens_avthetas = nc.Ensemble1_Average(thetas_list)
     ens_press = nc.Ensemble1_Average(press_list)
               
     #now get the perturbations
     #wvelthetaperts_list = []
     theta_pert_sq_list = []
     for i in range(len(wvels_list)):  #TODO: this should be more modular, see nchap_class                  
         thetapert_rough = np.subtract(thetas_list[i], ens_avthetas)
         thetapert = np.zeros_like(thetapert_rough)
         [znum, ynum, xnum] = wvels_list[i].shape
         for j in range(znum):#something like this is done in statistics.f90, staggered grid!
             if j == 0:
                 thetapert[j,:,:] = thetapert_rough[j,:,:]
             else:
                 thetapert[j,:,:] = 0.5*np.add(thetapert_rough[j,:,:], thetapert_rough[j-1,:,:])
         wvelpert = wvels_list[i]

         slice_lev = np.where(np.abs(height - hflux) < 26)[0][0]        
         wvelthetapert = np.multiply(wvelpert, thetapert)
         wvelperts_list.append(wvelpert[slice_lev, :, :])       
         thetaperts_list.append(thetapert[slice_lev, :, :])          

         [upwarm, downwarm, upcold, downcold]=nc.Flux_Quad_Wvels(wvelpert, thetapert) #TODO: expand clas Get_Vars.. to include this          

         upwarm_list.append(upwarm) 
         downwarm_list.append(downwarm)
         upcold_list.append(upcold)
         downcold_list.append(downcold)
         wvelthetaperts_list.append(wvelthetapert)     
         
     #and ensemble average them     
     #ens_upwarm = nc.Ensemble1_Average(upwarm_list)
     #ens_downwarm = nc.Ensemble1_Average(downwarm_list)
     #ens_upcold = nc.Ensemble1_Average(upcold_list)
     #ens_downcold = nc.Ensemble1_Average(downcold_list)
     #ens_avwvelthetaperts = nc.Ensemble1_Average(wvelthetaperts_list)
     #horizontally average them
     #upwarm_bar = nc.Horizontal_Average(ens_upwarm)     
     #downwarm_bar = nc.Horizontal_Average(ens_downwarm)
     #upcold_bar = nc.Horizontal_Average(ens_upcold)
     #downcold_bar = nc.Horizontal_Average(ens_downcold)
     #wvelthetapert_bar = nc.Horizontal_Average(ens_avwvelthetaperts)
     #applying new average, which ignores zero values
     upwarm_bar = nc.pert_h_Average(upwarm_list)     
     downwarm_bar = nc.pert_h_Average(downwarm_list)
     upcold_bar = nc.pert_h_Average(upcold_list)
     downcold_bar = nc.pert_h_Average(downcold_list)
     wvelthetapert_bar = nc.pert_h_Average(wvelthetaperts_list)               
     #save text files
     #print "SAVING", "/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/flux_quads_theta1" + dump_time 
     #np.savetxt("/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/flux_quads" + dump_time, np.transpose(np.array([upwarm_bar, downwarm_bar, upcold_bar, downcold_bar, wvelthetapert_bar])), delimiter=' ')
     #print upwarm_bar.shape, downwarm_bar.shape, upcold_bar.shape, downcold_bar.shape, wvelthetapert_bar.shape
     np.savetxt("/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/flux_quads_wvel1" + dump_time, np.transpose(np.array([upwarm_bar, downwarm_bar, upcold_bar, downcold_bar, wvelthetapert_bar])), delimiter=' ')
     
     #flatten the arrays, TODO: make a function or class method
     wvelperts = np.array(wvelperts_list)
     thetaperts = np.array(thetaperts_list)
     [enum, ynum, xnum] = wvelperts.shape

     #for a single case, and to look closer
     ##print 'where thetapert is less than -.5', np.where(thetaperts[0]<-.5)
     ##print 'wherre wvelpert is greater than .5', np.where(wvelperts[0]>.5)
     wvelperts_slice = wvelperts[0]
     thetaperts_slice = thetaperts[0]    

     wvelperts = np.reshape(wvelperts[0], ynum*xnum)
     thetaperts = np.reshape(thetaperts[0], ynum*xnum)
     wvelpertslice = wvelperts[0]
     thetapertslice = thetaperts[0]
     
     #wvelperts = np.reshape(wvelperts, enum*ynum*xnum)
     #thetaperts = np.reshape(thetaperts, enum*ynum*xnum)
     
     return height, wvelperts, thetaperts, wvelperts_slice, thetaperts_slice, upwarm_bar[slice_lev], downwarm_bar[slice_lev], upcold_bar[slice_lev], downcold_bar[slice_lev], wvelthetapert_bar[slice_lev]

go_ahead = np.int(input('have you changed the write out folder paths? 1 or 0: '))

if go_ahead == 1:
     
     date_list = ["Mar52014", "Jan152014_1", "Dec142013", "Nov302013", "Mar12014", "Dec202013", "Dec252013"] #"","", 
     
     #theFig2, theAxes2 = plt.subplots(nrows=3, ncols=3)     
     #theFig2.clf()
     #i=0
     #for theAx2 in theAxes2.flat:
     for i in range(len(date_list)):
         #print i
         #if i==2 or i==5:
         #    theAx2.axis('off')
         #else:         
         date = date_list[i]
         dump_time_list, Times = Make_Timelists(1, 1800, 28800)
         hvals = np.genfromtxt("/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/AvProfLims")
         #    scales = np.genfromtxt("/newtera/tera/phil/nchaparr/python/Plotting/"+date+"/data/invrinos")
         #    thetastar, wstar = scales[29, 9], scales[29, 2]
         #    lev_index = np.int(raw_input('which height level, 0, 1 or 2 (h0, h or h1)?:'))             
         #set up plots
         
         #theAx2 = theAxes2.flat[i]
             #theAx2.set_title(date, fontsize= 16) 
             #theAx2.set_title(r"$2d \ Histogram \ of \ Flux \ Quadrants$", fontsize= 16)
         for j in range(16):
             if date=="Nov302013":
                 hindex=2*(j+1)-1
             else:
                 hindex=3*(j+1)-1
         #    if i == 19:
         #need index specialized for Nov runs
             height, wvelperts, thetaperts, wvelperts_slice, thetaperts_slice, upwarm, downwarm, upcold, downcold, avflux = Main_Fun(date, dump_time_list[j], hvals[hindex, 1])
             ##print "Heights", hvals[29, 0], hvals[29, 1], hvals[29, 2], dump_time_list[19], Times[19]
               
         #av_quad_profs = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/flux_quads" + dump_time_list[i])

        #2d Hist
             #cmap = cm.hot
        #Estimate the 2D histogram
             #nbins = 200
             #H, xedges, yedges = np.histogram2d(1.0*wvelperts, 1.0*thetaperts, bins=nbins) #/wstar,/thetastar  
         #H needs to be rotated and flipped
             #H = np.rot90(H)
             #H = np.flipud(H)
        # Mask zeros
             #Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
         #theAx2.plot([-1, 1], [-1, 1])
        # Plot 2D histogram using pcolor
             #im = theAx2.pcolormesh(xedges,yedges,Hmasked, vmin = 0, vmax = 120, cmap = cmap)
                 #cbar = theFig2.colorbar(im)
                 #cbar.ax.set_ylabel(r'$Counts$')
             #theAx2.spines['left'].set_position('zero')
             #theAx2.spines['right'].set_color('none')
             #theAx2.spines['bottom'].set_position('zero')
             #theAx2.spines['top'].set_color('none')
             #theAx2.xaxis.set_ticks_position('bottom')
             #theAx2.yaxis.set_ticks_position('left')
             #theAx2.text(4, .5, r"$ w^{\prime} $ ",  fontdict=None, withdash=False, fontsize = 16)
             #theAx2.text(.5, 1, r"$ \theta^{\prime} $ ",  fontdict=None, withdash=False, fontsize = 16)
             #theAx2.set_ylim(-25, 25)
             #theAx2.set_xlim(-2, 3)
           
             #theAx2.set_ylim(-2, 1.5)
             #theAx2.set_xlim(-3, 5)
         i = i +1                    
            #theFig3.canvas.draw()
            #theFig1.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/fluxquadprofs.png")
            #theFig1.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/fluxquads.png")
            #theFig2.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/fluxquadhist"+str(lev_index)+".png")
            #theFig3.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/theta_cont"+str(lev_index)+".png")
     #theFig2.subplots_adjust(right=0.8)
     #cbar_ax = theFig2.add_axes([0.85, 0.15, 0.025, 0.7])
     #cbar_ax.set_ylabel(r'$Counts$')
     #theFig2.colorbar(im, cax=cbar_ax)
     #plt.show()

else:
    print('need to update write out folders') 
     


    
    
