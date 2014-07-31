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


ncFiles=glob.glob("/tera2/nchaparr/Dec142013/runs/sam_case1/OUT_3D/keep/NCHAPP1*_0000003600.nc") #3 D snaphsot files that have been converted to nc files
print ncFiles
"""
   For locating bug in Dec142013, ie getting  a key error for TABS

   Possible canditate for deletion 
"""
#set up plot
#theAx = nc.Do_Plot(1, 'Potential Temperature Profile Evolution', 'Height (m)', 'Potential Temperature (K)', 111)

index = 1
#loop over nc snapshot files
for theFile in ncFiles:
    print theFile
    ncdata = Dataset(theFile,'r')
    for (inname,invalue) in ncdata.variables.items():
    #create the variable
         print "printing varnames and values", inname, invalue
    
    #press = ncdata.variables['p'][...]
    #height = ncdata.variables['z'][...]
    #temp = np.squeeze(ncdata.variables['TABS'][...])
    #tracer = np.squeeze(ncdata.variables['TRACER'][...])
    ncdata.close()
    #meanpress = press
    #meantemp = np.mean(temp, axis = 1)#get the horizontal mean temperature
    #meantemp = np.mean(meantemp, axis = 1)
    #thetafact = np.array([(1.0*1000/i)**0.286 for i in meanpress]) #calculate potential temperature
    #theta = np.multiply(meantemp, thetafact)
    #print 'theta 0 from out3d', theta[0]
    #theAx.plot(theta, height) #, label = str(theFile)
    #plt.legend(loc = 'upper left', prop={'size':8})
    index = index + 1


#array1 = np.genfromtxt('initial1.txt') #get t0 dump from SAM
#array2 = np.genfromtxt('/tera/phil/nchaparr/python/Pert_Files/snd')
#initial_height2 = array2[:,0]
#initial_theta2 = array2[:,1]
#print 'theta 0 from snd', initial_theta2[0]
#initial_height = array1[:,0]
#initial_theta = array1[:,6]

#ncFile1= "/tera/phil/nchaparr/sam_ensemble/sam_case1/NCHAPP1/nchap1.nc" #initial nc file

#theAx.plot(initial_theta, initial_height,'*', label = 'Initial Sounding from setdata.f90')
#theAx.plot(initial_theta2, initial_height2,'bs', label = 'Initial Sounding from snd')

#plt.legend(loc = 'upper left', prop={'size':8})
#nc.Plot_nc(ncFile1, 1, theAx) #getting absolute temperatures from initial
                 #nc file and converting to potential temperature
#plt.legend(loc = 'upper left', prop={'size':8})
#plt.ylim(0, 2500)
#plt.xlim(295, 310)
#plt.show()



    
    
