import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from matplotlib import rcParams
rcParams.update({'font.size': 10})
import pandas as pd
import glob
from collections import defaultdict
import re
import os

#split 'theta_bar0000003420' into two parts
expr=re.compile('(.*)?(\d{10,10})')
file_dict=defaultdict(dict)
root_dir='/newtera/tera/phil/nchaparr/python/Plotting'
varnames=['theta_bar','press','wvelthetapert']
run_name_list = ["Nov302013","Dec142013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]
for case in run_name_list:
    for var_name in varnames:
        glob_expr='{}/{}/data/{}*'.format(root_dir,case,var_name)
        files=glob.glob(glob_expr)
        filenames={}
        for file_name in files:
            dirs,name=os.path.split(file_name)
            varname,the_time=expr.match(name).groups(1)
            filenames[int(the_time)]=file_name
        file_dict[case][var_name]=filenames
        print(case,var_name,len(file_dict[case][var_name]))

#all runs have the same heights, choose one at random
height_file='/newtera/tera/phil/nchaparr/python/Plotting/Mar52014/data/heights0000013200'
heights=np.genfromtxt(height_file)
num_levels,=heights.shape

time_full=np.linspace(600,28800,48)
time_nov302013=np.linspace(900,28800,32)


for case in run_name_list:
    for var_name in varnames:
        df=pd.DataFrame(heights,columns=['height'])
        if case == "Nov302013":
            time_list=time_nov302013
        else:
            time_list=time_full
        for the_time in time_list:
            df[the_time]=np.genfromtxt(file_dict[case][var_name][the_time])
            file_dict[case][var_name]['df']=df
#
# now do derived variables
#
for case in run_name_list:
    if case == "Nov302013":
        time_list=time_nov302013
    else:
        time_list=time_full
        for the_time in time_list:


with pd.HDFStore('paper_table.h5','r') as store:
     print(store.keys())
     df_cases=store.get('cases')

group_attributes={}
history='written 2015/8/5 by plot_theta_profs.py  9542a821e'        
group_attributes['/']=dict(time600=time_full,time900=time_nov302013,history=history,
                           case_list=run_name_list)

with pd.HDFStore('all_profiles.h5','w') as store:
    #loop over text files files
    store.put('df_overview',df_cases,format='table')
    #
    # write the vertical profiles
    # 
    for case in run_name_list:
        for var_name in varnames:
            node_name='/{}/{}'.format(case,var_name)
            store.put(node_name,file_dict[case][var_name]['df'],format='table')            

#     for i in range(len(theta_file_list)):
#         AvProfVars = np.genfromtxt(AvProfVars_list[i])
#         invrinosVars = np.genfromtxt(invrinos_list[i])
#         print(run_name,' invrinos ',invrinosVars.shape)
#         if AvProfVars.shape[0]==32:
#             big_array=np.empty([48,8],dtype=np.float64)
#             big_array[:,:]=np.nan
#             big_array[0:32,:]=AvProfVars[:,:]
#             AvProfVars=big_array[:,:]
#         if invrinosVars.shape[0]==32:
#             big_array=np.empty([48,10],dtype=np.float64)
#             big_array[:,:]=np.nan
#             big_array[0:32,:]=invrinosVars[:,:]
#             invrinosVars=big_array[:,:]
#         print('case: ',run_name,' shape: ',AvProfVars.shape)
#         columns=['h0','h','h1','zf0','zf','zf1','deltatheta','mltheta']
#         df_lims=pd.DataFrame(AvProfVars,columns=columns)
#         columns=['rino', 'invrino', 'wstar', 'S', 'tau', 'mltheta', 'deltatheta', 
#                  'pi3', 'pi4', 'thetastar']
#         df_rinos=pd.DataFrame(invrinosVars,columns=columns)
#        #Now for the gradients
#         dheight = np.diff(height)
#         dtheta = np.diff(theta)      
#         dthetadz = np.divide(dtheta, dheight)        
#         element0 = np.array([0])
#         dthetadz=np.hstack((element0, 1.0*dthetadz)) #*1.0/gamma
#         df_prof['dthetadz']=dthetadz
#         node_name='cases/{}/df_lims'.format(run_name)
#         store.put(node_name,df_lims,format='table')
#         node_name='cases/{}/df_prof'.format(run_name)
#         store.put(node_name,df_prof,format='table')
#         node_name='cases/{}/df_rinos'.format(run_name)
#         store.put(node_name,df_rinos,format='table')


#         #only need up to 2500meters
#         top_index = np.where(abs(1670 - height) < 40.)[0][0]

#         #where gradient is max, and flux is min
#         ##print AvProfVars[:,1].shape, height.shape

#         if run_name == "Nov302013":
#             h1 = AvProfVars[dump_time_index0, 1]
#         else:
#             h1 = AvProfVars[dump_time_index, 1]

#         h_index=np.where(dthetadz - np.amax(dthetadz[:top_index])==0)[0][0]
#         h=height[h_index]    
#         scaled_height = [1.0*ht/h for ht in height]

#         #print h1, h_index, height[h_index]

#         fluxes = np.multiply(wvelthetapert, rhow)*1004.0/flux_list[i]

#         Ax.plot(dthetadz, scaled_height, marker_list[i], label = legend_list[i], markersize=10) #, 

# zeros = np.zeros_like(height)
# Ax.plot(zeros+.01, scaled_height, 'k-')
# Ax.plot(zeros+.005, scaled_height, 'k-')
# Ax.plot(zeros+.005, scaled_height, 'k-')
# Ax.plot(zeros+.0025, scaled_height, 'k-')
# Ax.plot(zeros+.0002, scaled_height, 'k-')
# #Ax.legend(numpoints=1, loc = 'lower right', prop={'size':14})
# Ax.set_xticks([0.0002, .0025, .005, .010])
# Ax.set_xticklabels([0.2, 2.5, 5, 10])
# Ax.tick_params(axis="both", labelsize=25)
# plt.tight_layout()
# plt.show()
# #Fig1.savefig('/tera/phil/nchaparr/python/Plotting/Dec252013/pngs/theta_profs2hrs.png')
# #Fig1.savefig('/tera/phil/nchaparr/python/Plotting/Dec252013/pngs/flux_profs2hrs.png')





    
    
