import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import Make_Timelists
from nchap_class import For_Plots
import nchap_fun as nc
from matplotlib.lines import Line2D
import pdb

from matplotlib import rcParams
rcParams.update({'font.size': 10})

"""
   plots GM eq 26
  
"""

dump_time_list, Times = Make_Timelists(1, 600, 28800)
Times = np.array(Times)
dump_time_list0, Times0 = Make_Timelists(1, 900, 28800)
Times0 = np.array(Times0)

#plot


#Main Part -- pulling points and plotting them
label_list = ['100/10', '100/5', '60/5', '60/2.5', '150/5', '60/10', '150/10']
legend_list = ['kv', 'ko', 'yo', 'y*', 'ro', 'yv', 'rv']
Run_Date_List = ["Dec142013", "Nov302013", "Dec202013", "Dec252013", "Jan152014_1", "Mar12014", "Mar52014"]
keep_x = np.zeros([7,48],dtype=np.float)
keep_y = np.zeros_like(keep_x)

for i in range(len(label_list)):    
      
    if i<99:
        
        #pull variables from text files 
        points = For_Plots(Run_Date_List[i])
        rinovals = points.rinovals()
        gm_vars = points.gm_vars()
        heights = points.AvProfVars()
        
        #name and calculate points for plot
        z1_GM = heights[:,8]
        z_g=heights[:,1] # [elbot_dthetadz, h, eltop_dthetadz, elbot_flux, h_flux, eltop_flux, deltatheta, mltheta, z1_GM]
        delta = z1_GM - z_g
        z_enc=gm_vars[:, 3]
        L0 = gm_vars[:, 0]
        c_squared_delta = rinovals[:, 10]
        c_delta = c_squared_delta**.5

        Xs = (L0/z_enc) 
        Ys = delta/z_enc
        #Ys = points.Deltah_over_h(8, 1, 7)

        if Run_Date_List[i] == "Nov302013":
            #clean up
            #Xs=np.array(Xs)
            #Ys=np.array(Ys)
            #print(Xs, Xs.shape)
            #print(Ys, Ys.shape)
            #Ys[7:] = np.nan
            #Ys[13] = np.nan
            #Ys[15:17] = np.nan
            #Ys[24:26] = np.nan
            #
            # make all arrays the same size
            #
            Xs = np.hstack((Xs,np.empty(16,)))
            Ys = np.hstack((Ys,np.empty(16,)))
            Ys[31:48] = np.nan
            Xs[31:48] = np.nan
            #print(Xs, Xs.shape)
            #print(Ys, Ys.shape)
            #print(keep_x[i,:], keep_x[i,:].shape)
            #print(keep_y[i,:])
            #Xs[32:] = np.nan
            #Ys[32:] = np.nan
            #Ys[0] = np.nan #neg value
            #print('Nov 30 deltas: ',Ys)
            #plot 
        elif Run_Date_List[i] == "Jan152014_1":
             #clean up
            Ys[16:21] = np.nan
            Ys[:10] = np.nan
            Ys[29:] = np.nan
             #plot
        elif Run_Date_List[i] == "Mar12014":
            #clean up
            Ys[11:17] = np.nan
            Ys[:10] = np.nan
            Ys[29:] = np.nan
        else:
            Ys[:11] = np.nan
        
        keep_x[i,0:48] = Xs[0:48]
        keep_y[i,0:48] = Ys[0:48]
        #print(Run_Date_List[i])


plt.close('all')

fig1,ax1 = plt.subplots(1,1)
        
for count,name in enumerate(Run_Date_List):
    #if name=="Nov302013":
        #print(keep_x[count,:])
        #print(count)
        #print(keep_y[count,:])
        #print(keep_x[count,:])
    ax1.loglog(keep_x[count,:],keep_y[count,:], legend_list[count], label = label_list[count], markersize=10)

xs = np.arange(0, .2, .00001)
ys= (1)*xs**(.78)
ax1.loglog(xs, ys, 'k--')

ax1.text(.08, .15, r'$y = x^{.78}$',  fontdict=None, withdash=False, fontsize = 25, rotation=32)
ax1.set_xlabel(r"$(z_{enc}/L_{0})^{-1}$", fontsize=20)
ax1.set_ylabel(r"$\delta/z_{enc}$", fontsize=20)
plt.ylim(0.03, 0.2)
ax1.set_yticks([.1, .2])
ax1.set_yticklabels([0.1, .2])
plt.xlim(.03,.2)
ax1.set_xticks([.1, .2])
ax1.set_xticklabels([0.1, .2])
ax1.tick_params(axis="both", labelsize=20)
plt.tight_layout()


fig2,ax2 = plt.subplots(1,1)

x_keep=[]
y_keep=[]
for count,name in enumerate(Run_Date_List):
    #good_points= np.logical_not(np.isnan(keep_y[count,:]))
    #good_points=keep_x[count,:]
    good_x = keep_x[count,0:48]#[count,good_points]
    good_y = keep_y[count,0:48]#[count,good_points]
    hit = good_y < 0
    if np.sum(hit) > 0:
        print("found neg delta for {}".format(name))
    x=np.log10(good_x)
    y=np.log10(good_y)   
    ax2.plot(x,y, legend_list[count], label = label_list[count], markersize=10)
    x_keep.extend(good_x[0:48])
    y_keep.extend(good_y[0:48])

print(len(x_keep), len(y_keep))    
x_keep=np.array(x_keep)
y_keep=np.array(y_keep)
print(x_keep.shape, y_keep.shape)
fig2,ax2 = plt.subplots(1,1)

x=np.log10(x_keep)
y=np.log10(y_keep)
ax2.plot(x,y, legend_list[count], label = label_list[count], markersize=10)

from scipy.stats import linregress

logx=np.log10(x_keep)
logy=np.log10(y_keep)
mask=~np.isnan(logx)&~np.isnan(logy)
results = linregress(logx[mask], logy[mask])
xs = np.log10(np.arange(6, 40, 1))
ys=  results.intercept + results.slope*xs
ax2.plot(xs, ys, 'k--')
ax2.grid(True)

print('slope = ',results.slope)
print('intercept = ',10**results.intercept)

#pdb.set_trace()

plt.show()





    
    
