from netCDF4 import Dataset
import glob, os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import site
from datetime import datetime
#import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import nchap_fun as nc
from matplotlib.colors import Normalize
from Make_Timelist import *
from nchap_class import Get_Var_Arrays1
import sys
"""  
     For plotting 2d histograms of theta and wvel perturbations    
"""


def Main_Fun(date, dump_time, height_level):
    """Pulls output from an ensemble cases, gets ensemble averages and perturbations

    Arguments:
    date -- run date eg "Nov302013"
    dump_time -- time of output eg '0000000720'
    height_level -- height at which to take slice of perturbations, eg z_g0

    Returns:
    var_bar -- 64 array of horizontally averaged, ensemble averages or perturbations (covariances)
    
    """

    #create lists for variable arrays from each case
    wvelperts_list = []
    thetaperts_list = []

    #get velocity perts and thetas using method from nchap_class
    #sample file name for a case1
    #/tera/phil/nchaparr/tera2_cp/nchaparr/Aug122014/runs/sam_case1/OUT_3D/NCHAPP1_testing_doscamiopdata_24_0000000060.nc
    #
    filepath = {'root_dir': '/tera/phil/nchaparr/tera2_cp/nchaparr',
                'sam_dir': '/runs/sam_case'}
    suffix = "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_"
    filepath['date'] = date
    prefix = "{root_dir:}/{date:}/{sam_dir:}".format_map(filepath)
    Vars = Get_Var_Arrays1(prefix, suffix, dump_time)
    file_list = Vars.get_filelist()
    print(file_list)
    thetas_list, press_list = Vars.get_thetas()

    wvels_list = Vars.get_wvelperts()

    height = Vars.get_height()

    #get arrays of ensemble averaged variables, from nchap_fun

    ens_avthetas = nc.Ensemble1_Average(thetas_list)
    ens_press = nc.Ensemble1_Average(press_list)

    #now get the perturbations
    wvelthetaperts_list = []
    for i in range(len(
            wvels_list)):  #TODO: this should be more modular, see nchap_class
        thetapert_rough = np.subtract(thetas_list[i], ens_avthetas)
        thetapert = np.zeros_like(thetapert_rough)
        [znum, ynum, xnum] = wvels_list[i].shape
        for j in range(
                znum):  #something like this is done in statistics.f90, staggered grid!
            if j == 0:
                thetapert[j, :, :] = thetapert_rough[j, :, :]
            else:
                thetapert[j, :, :] = 0.5 * np.add(thetapert_rough[j, :, :],
                                                  thetapert_rough[j - 1, :, :])
        wvelpert = wvels_list[i]

        slice_lev = np.where(np.abs(height - height_level) < 26)[0][0]

        wvelthetapert = np.multiply(wvelpert, thetapert)
        wvelperts_list.append(wvelpert[slice_lev, :, :])
        thetaperts_list.append(thetapert[slice_lev, :, :])

        wvelthetaperts_list.append(wvelthetapert)

    #flatten the arrays, TODO: make a function or class method
    wvelperts = np.array(wvelperts_list)
    thetaperts = np.array(thetaperts_list)
    [enum, ynum, xnum] = wvelperts.shape

    #wvelperts_slice = wvelperts[0]
    #thetaperts_slice = thetaperts[0]

    wvelperts = np.reshape(wvelperts, enum * ynum * xnum)
    thetaperts = np.reshape(thetaperts, enum * ynum * xnum)

    return wvelperts, thetaperts


if __name__ == "__main__":
    # go_ahead = np.int(input('have you changed the read paths for hvals and scales? (yes=1 or no=0): '))
    # lev_index = np.int(input('which height level index in AvProfLims (z_f0=3, z_g0=0)?:'))

    go_ahead = 1
    lev_index = 0

    date_list = ["Mar52014", "Jan152014_1", "Dec142013", "Nov302013",
                 "Mar12014", "Dec202013", "Dec252013"]
    label_list = ["150/10", "150/5", "100/10", "100/5", "60/10", "60/5",
                  "60/2.5"]

    nov30_dump_time_list, nov30_Times = Make_Timelists(1, 900, 28800)
    nov30_time_index = 31

    dump_time_list, Times = Make_Timelists(1, 600, 28800)
    time_index = 47

    avproflims_prefix = "/tera/users/nchaparr/"
    avproflims_dir = "/data/AvProfLims"
    invrinos_dir = "/data/invrinos"

    the_dict = {}
    for date, label in zip(date_list, label_list):
        if date == "Nov302013":
            dump_time_list, Times = nov30_dump_time_list, nov30_Times
            time_index = nov30_time_index
        #get heights and convective scales from text files
        hval_file = "{}{}{}".format(avproflims_prefix, date, avproflims_dir)
        hvals = np.genfromtxt(hval_file)
        invrino_file = "{}{}{}".format(avproflims_prefix, date, invrinos_dir)
        scales = np.genfromtxt(invrino_file)
        thetastar, wstar = scales[time_index, 9], scales[time_index, 2]
        wvelperts, thetaperts = Main_Fun(date, dump_time_list[time_index],
                                         hvals[time_index, lev_index])
        the_dict[date] = (wvelperts, thetaperts, thetastar, wstar)

    #Get the perturbations using Main_Fun

    do_plots = False
    if do_plots:

        theFig2, theAxes2 = plt.subplots(nrows=3, ncols=3)

        #Loop over subplots for each date
        i = 0
        for theAx2 in theAxes2.flat:
            print('case', i)
            if i == 2 or i == 5:
                theAx2.axis('off')
            else:
                #Set up the axis spines
                theAx2.spines['left'].set_position('zero')
                theAx2.axvline(linewidth=2, color='k')
                theAx2.spines['right'].set_visible(False)
                theAx2.yaxis.set_ticks_position('left')
                theAx2.spines['bottom'].set_position('zero')
                theAx2.axhline(linewidth=2, color='k')
                theAx2.xaxis.set_ticks_position('bottom')
                theAx2.spines['top'].set_visible(False)

                #Set axis limits
                theAx2.set_yticks([-6, -3, 3, 6])
                theAx2.set_xticks([-2, -1, 1, 2])
                theAx2.tick_params(axis="both",
                                   length=10,
                                   width=2,
                                   direction='in')
                theAx2.tick_params(axis="x", pad=15)

                #Annotate
                theAx2.text(-2.4,
                            2,
                            r"$ \frac{w^{\prime}}{w^{*}} $ ",
                            fontdict=None,
                            withdash=False,
                            fontsize=16)
                theAx2.text(0.4,
                            -4,
                            label,
                            fontdict=None,
                            withdash=False,
                            fontsize=12)
                theAx2.text(.4,
                            4,
                            r"$ \frac{\theta^{\prime}}{\theta^{*}} $ ",
                            fontdict=None,
                            withdash=False,
                            fontsize=16)
                theAx2.set_ylim(-6, 6)
                theAx2.set_xlim(-3, 3)
                theAx2.set_axis_bgcolor('white')

                #Estimate the 2D histogram
                nbins = 100
                H, xedges, yedges = np.histogram2d(1.0 * wvelperts / wstar,
                                                   1.0 * thetaperts /
                                                   thetastar,
                                                   bins=nbins,
                                                   normed=True)

                #H needs to be rotated and flipped
                H = np.rot90(H)
                H = np.flipud(H)

                # Mask zeros
                Hmasked = np.ma.masked_where(
                    H == 0, H)  # Mask pixels with a value of zero

                # Plot 2D histogram
                xmin, xmax = np.amin(Hmasked), np.amax(Hmasked)
                vmin = 0
                vmax = .3
                the_norm = Normalize(vmin=vmin, vmax=vmax, clip=False)

                #Plot
                im = theAx2.pcolormesh(xedges,
                                       yedges,
                                       Hmasked,
                                       cmap='bone',
                                       alpha=0.5,
                                       norm=the_norm)  #, vmin = 0, vmax = 120

            #next subplot
            i = i + 1

        #Format Figure with colorbar
        theFig2.subplots_adjust(right=0.8)
        cbar_ax = theFig2.add_axes([0.85, 0.15, 0.025, 0.7])
        cbar_ax.set_xlabel(r"P($w^{\prime}, \theta^{\prime}$)/Bin Area",
                           fontsize=14)
        cbar_ax.xaxis.labelpad = 20
        cbar_ax.set_yticks([0, .1])
        cbar_ax.set_yticklabels([0, .1], fontsize=16)
        theFig2.colorbar(im, cax=cbar_ax)
        plt.show()
