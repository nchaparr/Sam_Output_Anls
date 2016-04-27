import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
rcParams.update({'font.size': 10})

"""

   Plots deltah_gamma for submission.

"""

date = "Mar52014"
sfc_flx = 150
gamma = .01

Fig1 = plt.figure(1)

#going to annotated using tex
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#two different sized subplots
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
Ax = Fig1.add_subplot(gs[0])
Ax1 = Fig1.add_subplot(gs[1])

Ax1.set_xlabel(r"$\overline{\theta}$", fontsize=33, labelpad=10)
Ax.set_xticks([])
Ax.set_xticklabels([])
Ax1.set_xticks([])
Ax1.set_xticklabels([])

Ax.set_ylabel(r"$z$", fontsize=33)
Ax.set_xlim(306.5, 310.5)
Ax1.set_xlim(306.5, 310.5)
dump_time_list, Times = Make_Timelists(1, 600, 28800)
 
theta_file_list = ["/tera/users/nchaparr/"+date+"/data/theta_bar"+ dump_time for dump_time in dump_time_list]
height_file = "/tera/users/nchaparr/"+date+"/data/heights0000000600"
AvProfVars = np.genfromtxt("/tera/users/nchaparr/"+date+"/data/AvProfLims")

i=37

theta = np.genfromtxt(theta_file_list[i])
height = np.genfromtxt(height_file)
        
#only need up to 2500meters
top_index = np.where(abs(2545 - height) < 40.)[0][0]

array = np.genfromtxt('/tera/users/nchaparr/Pert_Files/snd')
    
height_0 = array[:, 0]
theta_0 = array[:, 1]
f=interp1d(height_0, theta_0)

#Now plot inital sounding
top_index = np.where(height <= 2500)[0][-1]
theta_0 = f(height[0:top_index])

h0=AvProfVars[i,0]
h=AvProfVars[i,1]
h1=AvProfVars[i,2]
        
h0_index=np.where(height==h0)[0]
h_index=np.where(height==h)[0]
h1_index=np.where(height==h1)[0]
    
Ax.plot(theta, height, 'k-')
Ax1.plot(theta, height, 'k-', label = r"$\overline{\theta}$")#

Ax.annotate('', xy=(theta_0[h_index]-.1, h), xycoords = 'data', xytext=(theta[h1_index], h), textcoords = 'data', arrowprops=dict(arrowstyle = '<->'))
Ax.text(309.6, h-50, r"$\gamma \times \delta$", size=30)
Ax.annotate('', xy=(theta[h1_index], h), xycoords = 'data', xytext=(theta[h1_index], h1), textcoords = 'data', arrowprops=dict(arrowstyle = '<->'))
Ax.text(theta[h1_index]-.2, h+30, r"$\delta$", size=30)
Ax.set_yticks([h0, h, h1])
Ax.set_yticklabels([r"$z_{g0}$", r"$z_{g}$", r"$z_{g1}$"])
Ax1.set_yticks([25])
Ax1.set_yticklabels([0])
        
Ax.tick_params(axis="both", labelsize=30, width=3, length=15)
Ax1.tick_params(axis="both", labelsize=30, width=0, length=0)
Ax.set_ylim(h0-20, h1+20)
Ax1.set_ylim(25, 100)        

Ax.plot(theta_0 -.2, height[0:top_index], 'k--', label = r"$\overline{\theta}_{0}$") #,
Ax.spines['right'].set_visible(False)
Ax.spines['bottom'].set_visible(False)
Ax.spines['top'].set_visible(False)
Ax1.spines['right'].set_visible(False)
Ax1.spines['top'].set_visible(False)


d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=Ax.transAxes, color='k', clip_on=False)
Ax.plot((-d,+d),(-d,+d), **kwargs)# top-left diagonal
kwargs.update(transform=Ax1.transAxes)  # switch to the bottom axes
Ax1.plot((-d,+d),(1-3*d,1+3*d), **kwargs)   # bottom-left diagonal
print("about to show")
#plt.tight_layout()
plt.show()
#Fig1.savefig("/tera/phil/nchaparr/python/Plotting/"+date+"/pngs/theta_flux_profs.png")





    
    
