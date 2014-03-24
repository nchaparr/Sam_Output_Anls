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

def Main_Fun(daate, dump_time, hflux):
     
     """Pulls output from an ensemble cases, gets ensemble averages and perturbations and
     their horizontal averages

    Arguments:
    dump_time -- time of output eg '0000000720'

    Returns:
    var_bar -- 64 array of horizontally averaged, ensemble averages or perturbations (covariances)
    
    """
     #create list of filenames for given dump_time     
     ncfile_list = ["/tera2/nchaparr/"+date+"/runs/sam_case" + str(i+1) + "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_" + dump_time + ".nc" for i in range(10)]

     #create lists for variable arrays from each case
     upwarm_list = []
     downwarm_list = []
     upcold_list = []
     downcold_list = []
     wvelperts_list = []
     thetaperts_list = []

     #get velocity perts and thetas
     Vars =  Get_Var_Arrays1("/tera2/nchaparr/"+date+"/runs/sam_case", "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_", dump_time)         
     thetas_list, press_list = Vars.get_thetas()     
     wvels_list = Vars.get_wvelperts()          
     height = Vars.get_height()
     
     #get arrays of enseble averaged variables
     ens_avthetas = nc.Ensemble1_Average(thetas_list)
     ens_press = nc.Ensemble1_Average(press_list)
               
     #now get the perturbations
     wvelthetaperts_list = []
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

          [upwarm, downwarm, upcold, downcold]=nc.Flux_Quad(wvelpert, thetapert) #TODO: expand clas Get_Vars.. to include this          

          upwarm_list.append(upwarm) 
          downwarm_list.append(downwarm)
          upcold_list.append(upcold)
          downcold_list.append(downcold)
          wvelthetaperts_list.append(wvelthetapert)     
          
     #and ensemble average them     
     ens_upwarm = nc.Ensemble1_Average(upwarm_list)
     ens_downwarm = nc.Ensemble1_Average(downwarm_list)
     ens_upcold = nc.Ensemble1_Average(upcold_list)
     ens_downcold = nc.Ensemble1_Average(downcold_list)         
     ens_avwvelthetaperts = nc.Ensemble1_Average(wvelthetaperts_list)

     #horizontally average them
     upwarm_bar = nc.Horizontal_Average(ens_upwarm)     
     downwarm_bar = nc.Horizontal_Average(ens_downwarm)
     upcold_bar = nc.Horizontal_Average(ens_upcold)
     downcold_bar = nc.Horizontal_Average(ens_downcold)
     wvelthetapert_bar = nc.Horizontal_Average(ens_avwvelthetaperts)
          
     #save text files
     np.savetxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/flux_quads" + dump_time, np.transpose(np.array([upwarm_bar, downwarm_bar, upcold_bar, downcold_bar, wvelthetapert_bar])), delimiter=' ')
     
     #flatten the arrays, TODO: make a function or class method
     wvelperts = np.array(wvelperts_list)
     thetaperts = np.array(thetaperts_list)
     [enum, ynum, xnum] = wvelperts.shape

     #for a single case, and to look closer
     #print 'where thetapert is less than -.5', np.where(thetaperts[0]<-.5)
     #print 'wherre wvelpert is greater than .5', np.where(wvelperts[0]>.5)
     wvelperts_slice = wvelperts[0]
     thetaperts_slice = thetaperts[0]    

     #wvelperts = np.reshape(wvelperts[0], ynum*xnum)
     #thetaperts = np.reshape(thetaperts[0], ynum*xnum)
     #wvelpertslice = wvelperts[0]
     #thetapertslice = thetaperts[0]
     
     wvelperts = np.reshape(wvelperts, enum*ynum*xnum)
     thetaperts = np.reshape(thetaperts, enum*ynum*xnum)
     
     return height, wvelperts, thetaperts, wvelperts_slice, thetaperts_slice, upwarm_bar[slice_lev], downwarm_bar[slice_lev], upcold_bar[slice_lev], downcold_bar[slice_lev], wvelthetapert_bar[slice_lev]

go_ahead = np.int(raw_input('have you changed the write out folder paths? 1 or 0: '))
if go_ahead == 1:

     date = "Dec252013"
     
     dump_time_list, Times = Make_Timelists(1, 600, 28800)
     hvals = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/AvProfLims")
     lev_index=2
        
     #set up plots
     #theFig = plt.figure(3)     
     #theFig.clf()
     #theAx = theFig.add_subplot(111)
     #theAx.set_title(r"$Flux \ Quadrants$", fontsize= 16)
     #theAx = nc.Do_Plot(3, r"$Flux \ Quadrants$", fontsize= 16, '', '', 111)
     #Todo: add option to take args to Do_Plot

     #theFig1 = plt.figure(4)     
     #theFig1.clf()
     #theAx1 = theFig1.add_subplot(111)
     #theAx1.set_title(r"$Flux \ Profiles$", fontsize= 16)
     #theAx1 = nc.Do_Plot(fignum, title, ylabel, xlabel, sub)

     #theFig2 = plt.figure(5)     
     #theFig2.clf()
     #theAx2 = theFig2.add_subplot(111)
     #theAx2.set_title(r"$2d \ Histogram \ of \ Flux \ Quadrants$", fontsize= 16)
     #theAx2 = nc.Do_Plot(fignum, title, ylabel, xlabel, sub)

     #for single case contours of theta, w
     theFig3 = plt.figure(2)     
     theFig3.clf()
     theAx3 = theFig3.add_subplot(111)
     #theAx3.set_title(r"$Contour \ of \theta^{,}$", fontsize= 16)
     
     #get horizontally averaged ensemble averaged variable and plot
     colorlist=['k', 'b', 'c', 'g', 'r', 'm', 'y', '.75']
     for i in range(48):
          if i == 11:
               
               height, wvelperts, thetaperts, wvelperts_slice, thetaperts_slice, upwarm, downwarm, upcold, downcold, avflux = Main_Fun(date, dump_time_list[i], hvals[i, lev_index])
               
               av_quad_profs = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/flux_quads" + dump_time_list[i])

               #theAx1.plot(av_quad_profs[:, 0], height,'r-', label = 'up warm')
               #theAx1.plot(av_quad_profs[:, 1], height, 'b--', label = 'down warm')
               #theAx1.plot(av_quad_profs[:, 2], height, 'b-', label = 'up cold')
               #theAx1.plot(av_quad_profs[:, 3], height, 'r--', label = 'down cold')
               #theAx1.plot(av_quad_profs[:, 4], height, 'k-', label = 'average')
               #theAx1.plot(np.zeros_like(height), height, 'k-')
               #theAx1.set_ylim(100, 2000)
               #theAx1.legend(loc = 'upper right', prop={'size':8})
               
               #theAx.plot(wvelperts, thetaperts, 'ro', markersize=1, markeredgecolor='none')
               #theAx.spines['left'].set_position('zero')
               #theAx.spines['right'].set_color('none')
               #theAx.spines['bottom'].set_position('zero')
               #theAx.spines['top'].set_color('none')
               
               #theAx.xaxis.set_ticks_position('bottom')
               #theAx.yaxis.set_ticks_position('left')
               #theAx.set_ylim(-1.5, 1.5)
               #theAx.set_xlim(-3, 5)
                              
               #theAx.text(2, 1, "$%.5f$"%upwarm,  fontdict=None, withdash=False, fontsize = 16)
               #theAx.text(2, -1, "$%.5f$"%upcold,  fontdict=None, withdash=False, fontsize = 16)
               #theAx.text(-2, 1, "$%.5f$"%downwarm,  fontdict=None, withdash=False, fontsize = 16)
               #theAx.text(-2, -1, "$%.5f$"%downcold,  fontdict=None, withdash=False, fontsize = 16)
               #theAx.text(3.5, .2, r"$ w^{,} $ ",  fontdict=None, withdash=False, fontsize = 16)
               #theAx.text(-.5, 1.25, r"$ \theta^{,} $ ",  fontdict=None, withdash=False, fontsize = 16)

               #2d Hist
               #cmap = cm.autumn
               # Estimate the 2D histogram
               #nbins = 200
               #H, xedges, yedges = np.histogram2d(wvelperts, thetaperts, bins=nbins)
                # H needs to be rotated and flipped
               #H = np.rot90(H)
               #H = np.flipud(H)
               # Mask zeros
               #Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
                # Plot 2D histogram using pcolor
               #im = theAx2.pcolormesh(xedges,yedges,Hmasked, vmin = 0, vmax = 120, cmap =cmap)
               #cbar = theFig2.colorbar(im)
               #cbar.ax.set_ylabel(r'$Counts$')
                                      
               #theAx2.spines['left'].set_position('zero')
               #theAx2.spines['right'].set_color('none')
               #theAx2.spines['bottom'].set_position('zero')
               #theAx2.spines['top'].set_color('none')
               #theAx2.xaxis.set_ticks_position('bottom')
               #theAx2.yaxis.set_ticks_position('left')
               #theAx2.text(3.5, .2, r"$ w^{,} $ ",  fontdict=None, withdash=False, fontsize = 16)
               #theAx2.text(-.5, 1.25, r"$ \theta^{,} $ ",  fontdict=None, withdash=False, fontsize = 16)
               
               #theAx2.set_ylim(-1.5, 1.5)
               #theAx2.set_xlim(-3, 5)

               #theAx3.set_title(r"$Contour \ of \ \theta^{,} \ after \ " + str(Times[i]) +"\ hours$")
               theAx3.set_xlabel(r"$x \ (m)$")
               theAx3.set_ylabel(r"$y \ (m)$")

               v_max, v_min, mean, stddev = np.amax(thetaperts_slice), np.amin(thetaperts_slice), np.mean(thetaperts_slice), np.std(thetaperts_slice)

               filler_array = np.zeros([64, 192])
               
               Slice = np.vstack((thetaperts_slice, filler_array))
               x = np.arange(0, 4800, 25)
               y = np.arange(0, 4800, 25)
               X,Y = np.meshgrid(x, y)
         
               im = theAx3.pcolor(X, Y, np.transpose(Slice), cmap=cm.hot, vmax=v_max, vmin=v_min)
               bar = theFig3.colorbar(im)
               theAx3.set_xlim(0, 3200)
               theAx3.set_ylim(0, 4800)
                              
               theFig3.canvas.draw()
               #theFig.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/fluxquadprofs0.png")
               #theFig1.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/fluxquads0.png")
               #theFig2.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/fluxquadhist0.png")
               theFig3.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/theta_cont"+str(lev_index)+".png")
     plt.show()

else:
    print 'need to update write out folders' 
     


    
    
