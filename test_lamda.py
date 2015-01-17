import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt

"""testing function get_level for use in get_rino and get_dheta
"""
F0 = 75 #j/sm2
rho = 1 #Kg/m3
cp = 1004 #j/KgK

#TODO: try generalizing what's done further down with this fn
#TODO: see if an inequality can, or if statement can be passed via a handel or lamda fn
def get_level(array,  condition1, val11, condition2, val12):
    """
    Arguments:
    array -- 
    val11, val12 -- 
    condition1, 2 -- lamda functions
    
    Returns:
    level -- index where condition is met  
    
    """
    
    for k in range(len(array)-1):            
         if condition1(array, k, val11) and condition2(array, k, val12):
              level = k+1
              print array[k+1], k+1
              break    
    return level


#create lists of txt file to loop over
Times  = [180, 360, 540, 720, 900, 1080]
dump_time_list = ['0180', '0360', '0540', '0720', '0900', '1080']
theta_file_list = ["/tera/phil/nchaparr/SAM2/sam_main/python/Plotting/May172013/theta_bar000000"+ dump_time for dump_time in dump_time_list]
flux_file_list = ["/tera/phil/nchaparr/SAM2/sam_main/python/Plotting/May172013/wvelthetapert000000"+ dump_time for dump_time in dump_time_list]
height_file = "/tera/phil/nchaparr/SAM2/sam_main/python/Plotting/May172013/heights0000000180"
press_time = np.genfromtxt('press')

for i in range(len(theta_file_list)):
    theta = np.genfromtxt(theta_file_list[i])
    height = np.genfromtxt(height_file)
    press = press_time[:,i]
    dheight = np.diff(height)
    dtheta = np.diff(theta)
    
    dthetadz = np.divide(dtheta, dheight)
    
    element0 = np.array([0])
    dthetadz=np.hstack((element0, dthetadz))
        
    #only need up to 2500meters
    top_index = np.where(abs(2500 - height) < 26.)[0][0]

    #where gradient is greater than zero    
    xgty = lambda array, i, y : array[i+1] > y
    xgtoey = lambda array, i, y :array[i] >= y    
    dtheta_index_b = get_level(dthetadz, xgty, .0002, xgtoey, 0.00)
             
    #where gradient resumes as gamma   
    xety = lambda array, i, y : abs(array[i+1]-y)<.0002
    xgty1 = lambda array, i, y: array[i] > array[i+y]
    dtheta_index_t = get_level(dthetadz, xety, .005, xgty1, 1)
            
    #now fluxes
    fluxes = np.genfromtxt(flux_file_list[i])

    xltoey = lambda array, i, y: array[i+1] <= y
    xgty2 = lambda array, i, y: array[i] > y
    flux_index_b = get_level(fluxes, xltoey, 0, xgty2, 0)
    
    xlty = lambda array, i, y: abs(array[i+1]) < y
    xlty1 = lambda array, i, y: array[i] < y
    flux_index_b = get_level(fluxes, xlty, 0.3, xlty1, 0)
    
    
    



    
    
