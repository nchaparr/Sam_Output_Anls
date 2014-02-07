from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
import matplotlib.pyplot as plt
import site
site.addsitedir('/tera/phil/nchaparr/python')
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from Make_Timelist import *
from tempfile import TemporaryFile


"""
   Gets horizontal slices at z levels
   for doing FFTs
   don't know if I need it
"""
dump_time_list, Times = Make_Timelists(1, 600, 28800)
print dump_time_list[20]
AvProfVars = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/Dec142013/data/AvProfLims")

for dump_time in dump_time_list:
     if dump_time == dump_time_list[20]:    
    
         [wvels, theta, tracer, height] = nc.Get_Var_Arrays('/tera2/nchaparr/Dec142013/runs/sam_case', '/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_', dump_time, 1)

         height_index = np.where(height - AvProfVars[20, 1] == 0)[0][0]
         print height_index

         wvelslice_h = wvels[height_index, :, :] 

         outfile = TemporaryFile()

         np.savez('/tera/phil/nchaparr/python/Plotting/Dec142013/data/wvelslice_h'+ dump_time + '.npz', wvelslice_h)

    
    
    



    
    
