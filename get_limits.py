from __future__ import print_function
from __future__ import absolute_import
import ruamel.yaml
import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Sam_Output_Anls.Make_Timelist import Make_Timelists
import sys
from Sam_Output_Anls import nchap_fun as nc
from matplotlib import rcParams
import os
from collections import OrderedDict as od
rcParams.update({'font.size': 10})



"""calculates temperature gradients (discrete) from txt files inturn from ensemble run 3D files
   gets levels where gradient exceeds zero, and where it resumes gamma, and zero crossings for fluxes
   and the maximal points withing the entrainment region.
   Gets Delta Theta and the Mixed Layer average Theta.
   Dumps them in a text file.
   Calcs and dumps rino, invrino, wstar
"""
#to be changed for each run

cases = """
   Nov302013,.005,100,39; Dec142013,   .01,100,23; Dec202013,.005,60,31;
   Dec252013,.0025,60,51; Jan152014_1,.005,150,48; Mar12014,  .01,60,18;
   Mar52014,  .01,150,29 
"""

write_yaml = True
if write_yaml:
    def sort_L0(item):
        return item['L0']

    ind_cases = cases.split(';')
    case_list = []
    key_list = ['name','gamma','flux','L0']
    for index,case in enumerate(ind_cases):
        case = case.strip()
        case_tup = case.split(',')
        val_dict=od()
        for key,value in zip(key_list,case_tup):
            if key != 'name':
                value = float(value)
            val_dict[key] = value
        if val_dict['name'] == 'Nov302013':
            val_dict['del_t'] = 900
        else:
            val_dict['del_t'] = 600
        case_list.append(val_dict)
    case_list.sort(key = sort_L0)

    case_file = 'cases.yaml'
    with open(case_file,'w') as f:
        ruamel.yaml.dump(case_list,f,Dumper=ruamel.yaml.RoundTripDumper,default_flow_style=False)
else:
    case_list = ruamel.yaml.load(f)
    

            
    
for index,case in enumerate(case_list):
    print('case number: ',index)
    rundate = case['name']
    gamma = case['gamma']
    flux_s = case['flux']

    #output times
    del_t = case['del_t']
    dump_time_list, Times = Make_Timelists(1, del_t, 28800)
    Times = np.array(Times)  

    input_dir = "/tera/users/nchaparr/"
    output_dir = './fig6c/'

    theta_file_list = [input_dir + rundate + "/data/theta_bar"+ dump_time for dump_time in dump_time_list]
    press_file_list = [input_dir + rundate + "/data/press"+ dump_time for dump_time in dump_time_list]
    flux_file_list = [input_dir + rundate + "/data/wvelthetapert"+ dump_time for dump_time in dump_time_list]
    height_file = input_dir + rundate + "/data/heights0000000600"

    AvProfLims = []
    invrinos = []
    #loop over text files files
    for i in range(len(theta_file_list)):

        theta = np.genfromtxt(theta_file_list[i])
        height = np.genfromtxt(height_file)

        press = np.genfromtxt(press_file_list[i])
        rhow = nc.calc_rhow(press, height, theta[0])
        wvelthetapert = np.genfromtxt(flux_file_list[i])

        #Now for the gradients
        dheight = np.diff(height)
        dtheta = np.diff(theta)

        dthetadz = np.divide(dtheta, dheight)

        element0 = np.array([0])
        dthetadz=np.hstack((element0, dthetadz))

        #only need up to 1900meters
        top_index = np.where(abs(1700 - height) < 26.)[0][0] #may need to be higher (e.g. for 60/2.5)

        #TODO: see test_lamda
        #where gradient is greater than zero    
        for j in range(len(dthetadz)-1):
            if (dthetadz[j+1] >.0002) and (dthetadz[j] >= 0):
                dtheta_index_b = j+1
                break

        #where gradient resumes as gamma    
        for k in range(len(dthetadz[:top_index])-1):        
            if np.abs(dthetadz[k+2]-gamma)<.0002 and np.abs(dthetadz[k+1]-gamma)<.0002 and dthetadz[k-1]>gamma:            
                dtheta_index_t = k+1                        
                break

        #now fluxes    
        fluxes = np.multiply(wvelthetapert, rhow)*1004.0
        for l in range(len(dthetadz)-1):
            if (fluxes[l+1] <= .0) and (fluxes[l] > 0):
                flux_index_b = l+1
                break

        for m in range(len(dthetadz[0:top_index])-1):
            if (abs(fluxes[m+1]) < 0.5) and (fluxes[m] < 0) and (fluxes[m-1] < 0):
                flux_index_t = m+1
                break

        mltheta = np.mean(theta[0:dtheta_index_b])
        deltatheta = -theta[dtheta_index_b]+theta[dtheta_index_t]
        eltop_dthetadz = height[dtheta_index_t]
        elbot_dthetadz = height[dtheta_index_b]
        eltop_flux = height[flux_index_t]
        elbot_flux = height[flux_index_b]
        h = height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]]
        h_flux = height[np.where(wvelthetapert - np.amin(wvelthetapert) == 0)[0][0]]

        #TODO: this can be tidied up, ie name valriables and pass the named variables to calc_rino    
        print(i, height[dtheta_index_b], height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]], height[dtheta_index_t], height[flux_index_b], height[np.where(wvelthetapert - np.amin(wvelthetapert) == 0)[0][0]], height[flux_index_t])
        print(i, height[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]], np.mean(theta[0:dtheta_index_b]), 1.0*flux_s/(rhow[0]*1004), -theta[dtheta_index_b]+theta[dtheta_index_t])

        [rino, invrino, wstar, S] =  nc.calc_rino(h, mltheta, 1.0*flux_s/(rhow[0]*1004), deltatheta, gamma)

        AvProfLims.append([elbot_dthetadz, h, eltop_dthetadz, elbot_flux, h_flux, eltop_flux, deltatheta, mltheta])
        tau = 1.0*h/wstar
        invrinos.append([rino, invrino, wstar, S, tau, mltheta, deltatheta])


    files = ['AvProfLims_old','invrinos_old']
    data_arrays = [np.array(AvProfLims),np.array(invrinos)]
    for the_file,the_array in zip(files,data_arrays):
        path = output_dir + rundate + '/data/'
        if not os.path.exists(path):
            os.makedirs(path)
        full_path = path + the_file
        np.savetxt(full_path,the_array, delimiter=' ') #' + rundate + '/
