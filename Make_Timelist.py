from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib
import matplotlib.pyplot as plt
import site
from . import nchap_fun

"""
   Get list of strings representing times of Output
   
"""

#TDOD: put into nchap_fun?

def make_string(num):
     """Integer to string

    Arguments:
     integer -- 

    Returns:
     string -- eg '0000000030'     
    
    """
     
     if num<100:
          string = '00000000' + str(num)

     if num >99 and num<1000:
          string = '0000000' + str(num)

     if num >999 and num<10000:
          string = '000000' + str(num)

     if num >9999 and num<100000:
          string = '00000' + str(num)

     if num >99999 and num<1000000:
          string = '0000' + str(num)     

     return string

def Make_Timelists(dt, time_diff, stop_time):
     """Main Function

    Arguments:
     dt -- from prm file
     time_diff -- interval between output writes
     stop_time -- end of run from prm
     
    Returns:
     dump_time_list, Times_hrs -- list of strings for filenames. list of times      
    
    """
     num_times = int(1.0*stop_time/time_diff)
          
     Times = [(i+1)*time_diff for i in range(num_times)]
     
     Times_hrs = [1.0*dt*Time/3600 for Time in Times]
          
     dump_time_list = [make_string(Time) for Time in Times]

     
     return dump_time_list, Times_hrs

if __name__ == "__main__":
     dump_time_list, Times_hrs = Make_Timelists(2, 30, 14400)
