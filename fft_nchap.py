from __future__ import division
import glob,os.path
import numpy as np
import scipy.fftpack as fftpack 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import site
from datetime import datetime
from nchap_class import *

site.addsitedir('/tera/phil/nchaparr/python')
#from nchap_fun import *
site.addsitedir('/tera/phil/nchaparr/python/Plotting/Dec142013/python')
from Make_Timelist import *
#import sys
#sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc


"""
   Starting to build script to plot power spectrum at a slice

   First mirroring what Phil did, except for the sf window
   but I will make a square slice     

"""

date="Dec252013"
#get the square slice from the nc file
dump_time_list, Times = Make_Timelists(1, 600, 28800)
ncfile_list = ["/tera2/nchaparr/"+date+"/runs/sam_case" + str(i+1) + "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_" + dump_time_list[11] + ".nc" for i in range(10)]
thefile = ncfile_list[0]
#thefile = "a17.nc"
ncdata = Dataset(thefile,'r')
#tau = ncdata.variables['tau'][...]
#print tau.shape
#[ynum, xnum] = tau.shape

data = For_Plots(date)
AvProfVars = data.AvProfVars()
RinoVals = data.RinoVals()
#AvProfVars = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/AvProfLims")
#RinoVals = np.genfromtxt("/tera/phil/nchaparr/python/Plotting/"+date+"/data/invrinos")

h0 = AvProfVars[11, 0]
h = AvProfVars[11, 1]
h1 = AvProfVars[11, 2]
wstar = RinoVals[11, 2]
height = np.squeeze(ncdata.variables['z'][...])

lev0 = np.where(abs(h0 - height) < 6.)[0][0]
#print lev0
wvel0 = np.squeeze(ncdata.variables['W'][...])[lev0,:,:]
[ynum, xnum] = wvel0.shape
wvel_slice0 = wvel0[:, :ynum]

lev1 = np.where(abs(h - height) < 6.)[0][0]
#print lev1
wvel1 = np.squeeze(ncdata.variables['W'][...])[lev1,:,:]
wvel_slice1 = wvel1[:, :ynum]

lev2 = np.where(abs(h1 - height) < 6.)[0][0]
#print lev2
wvel2 = np.squeeze(ncdata.variables['W'][...])[lev2,:,:]
wvel_slice2 = wvel2[:, :ynum]

mid_point = int(np.floor(1.0*ynum/2))
print mid_point
  

def anav(wvel_slice):
     #get the fft
     slice_fft = fftpack.fft2(wvel_slice)
     
     #shift fft
     
     slice_fft = fftpack.fftshift(slice_fft)
     
     #get the energy density function
     e_dense = slice_fft*np.conj(slice_fft)/(ynum*ynum)
     e_dense = e_dense.real
        
    #Anulus Averaging
     nbns = int(round((np.sqrt(2)*ynum/1),0)+1) #sqrt(ynum**2 + ynum**2)
     e_spec = np.zeros(nbns)
     cnt = np.zeros(nbns)
     for i in range(ynum):
         for j in range(ynum):
         
             r = np.sqrt(((i+1)-ynum/2)**2 + ((j+1)-ynum/2)**2)
       #      print i, j, r, np.sqrt(((i+1))**2 + ((j+1))**2)
             bn = int(np.floor(r/1))
      #       print bn
             e_spec[bn]=e_spec[bn]+r*e_dense[i,j]
             cnt[bn]=cnt[bn]+1
     for i in range(nbns):
         if cnt[i]>0:
             #e_spec[i]=e_spec[i]/cnt[i]/(4*(np.pi**2)) #why
             e_spec[i] = np.pi*e_spec[i]/cnt[i]
     #print nbns, mid_point
     
     x = np.arange(nbns)+1
     
     x_spec = x[0:mid_point] #get k in 1/meters? ie DeltaX, Deltay = 25 meters     
     delta_k = 1./.025                # 1./km 
     nyquist = delta_k*0.5
     num_delta_ks = len(x_spec)
     #print .025*num_delta_ks
     x_spec = x_spec*(delta_k/float(num_delta_ks))# k = w/(25m)
     
     e_spec = e_spec[0:mid_point]

     slope = x_spec**(-5/3) #why this not, -5/3?
     return e_dense, x_spec, e_spec, slope

e_dense0, x_spec0, e_spec0, slope0 = anav(wvel_slice0)
e_dense1, x_spec1, e_spec1, slope1 = anav(wvel_slice1)
e_dense2, x_spec2, e_spec2, slope2 = anav(wvel_slice2)


Fig = plt.figure(1)
Ax = Fig.add_subplot(111)
#Ax.set_title('')
Ax.loglog(x_spec0, e_spec0, 'ko', label = r'$h_{0}$') 
Ax.loglog(x_spec1, e_spec1, 'kv', label = r'$h$')
Ax.loglog(x_spec2, e_spec2, 'k*', label = r'$h_{1}$')
Ax.loglog(x_spec0[10:], slope0[10:]*5000,'r--', label='slope = -5/3')
Ax.loglog(x_spec0[10:], slope0[10:]*300,'r--')
Ax.legend(loc = 'upper left', prop={'size': 10}, numpoints=1)
#plt.xlim(.0, 10)
#Ax.set_ylim(10**-6, 10**2)
Ax.set_xlabel(r"$k / (1.6km)$", fontsize=20)
Ax.set_ylabel(r"$E(k)$", fontsize=20)
plt.show()

theFig3 = plt.figure(2)     
theFig3.clf()
theAx3 = theFig3.add_subplot(111)
theAx3.set_title(r"$E(k_{x}, k_{y})$", fontsize= 16)
theAx3.set_xlabel(r"$k_{x} \ (1/m)$")
theAx3.set_ylabel(r"$k_{y} \ (1/m)$")
freqs = fftpack.fftfreq(128, 1.0/128)
print freqs
k_array = np.arange(1, 65, 1)
#print k_array
k_array=k_array[::-1]
#print k_array
k_array = np.hstack((k_array, -np.arange(0, 64, 1)))
#print k_array, k_array.shape
X, Y = np.meshgrid(k_array, k_array)
     
e_dense0 = np.ma.masked_where(e_dense0==0, e_dense0)
v_max, v_min, mean, stddev = np.amax(e_dense0), np.amin(e_dense0), np.mean(e_dense0), np.std(e_dense0)
norm = matplotlib.colors.Normalize(vmin = np.amin(e_dense0), vmax = np.amax(e_dense0), clip = False)
im = theAx3.pcolor(X, Y, np.transpose(e_dense0), cmap=cm.hot, norm=norm)
bar = theFig3.colorbar(im)
#CS = plt.contour(X, Y, np.transpose(e_dense_masked))
theAx3.set_xlim(-32, 32)
theAx3.set_ylim(-32, 32)
plt.show  

#Fig.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/2dffts_v.png")
    
