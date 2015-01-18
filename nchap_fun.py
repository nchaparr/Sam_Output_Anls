from netCDF4 import Dataset
import glob,os.path
import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt

""" starting to collect commonly used functions"""

def calc_rino(BLHeight, MLTheta, SfcFlux, Theta_jump, gamma, delta_h):
     """Richardson number

    Arguments:
    BLHeight -- height of mixed layer, scalar, m
    MLTheta -- average potential temperature in mixed layer, scalar, K
    SfcFlx -- surface heat flux, scalar, w/m2
    Theta_jump -- inversion temperature jump, scalar, K

    Returns:
    rino, invrino, wstar -- Richardson Number, Inverse Richardson Number, convective velocity scale
    
    """
     intermed = 1.0*BLHeight/MLTheta
     intermed1 = intermed*SfcFlux
     wstar = (9.81*intermed1)**(1.0/3)
     thetastar = 1.0*SfcFlux/wstar
     rino = 1.0*Theta_jump/thetastar
     
     S = ((1.0*BLHeight/wstar)**2)*(gamma)*(1.0*9.81/MLTheta)

     pi3 = gamma*1.0*BLHeight/Theta_jump

     pi4 = gamma*1.0*delta_h/Theta_jump
     
     return rino, 1.0/rino, wstar, S, pi3, pi4


def get_dhdt(heights, time):
     #TODO: see if this an fn for calc ing dthetadz can be unified
     """

    Arguments:
    heights  -- array 
    time -- array

    Returns:
    dhdt -- array     
    
    """
     dh = np.diff(heights)
     dt = np.diff(time)*3600
     dhdt = np.divide(dh, dt)     
     return dhdt

def calc_rhow(press, height, ts):
     """hydrostatic relation

    Arguments:
    press, height -- arrays 

    Returns:
    ts -- surface temperature (K)     
    
    """
    
     rhow = np.zeros_like(press)
     

     rhow[0] = 1.0*press[0]*100/(287*ts)
     #TODO: check if this is ok
     
     for i in range(press.shape[0]-1):
         
         rhow[i+1] = 1.0*(press[i] - press[i+1])/(height[i+1] - height[i])*(1.0*100/9.81)
         
     return rhow

def calc_dh(DTheta, F0, rho, cp, dt):
    """Gets change in height for a change in time per phil's notes?

    Arguments:
    DTheta -- change in potential temperature, K
    F0 -- surface heat flux, w/m2
    rho -- density of air, kg/m3
    cp -- heat capacity of air
    dt -- change in time, s

    Returns:
    dh -- change in height, m  
    
    """     
    dh = 1.0*(dt)*((0.2*F0)/(rho*cp*DTheta))
    return dh

def abs2pot(press0, press, abstemp, reverse):#see SAM6.8.2/SRC/set.data.f90
    """Calculates array of potential temperatures from heights and absolute temperatures

    Arguments:
    press -- 1d array of pressures
    abstemp -- 1d array of absolute temperatures of height.shape
    reverse? -- 0 or 1, instruction to reverse order of arrays, ie for nc files

    Returns:
    pottemp --  1d array of potential temperatures 

    """
    end=len(press) #index for last element
        
    ggr =  9.81#m/s2
    cp = 1004 #j/Kg/K
    rgas = 287#J/kg/K
    
    height = []
    pottemp = []

    #reverse arrays for looping
    if reverse == 1:
        press = press[::-1]
        abstemp = abstemp[::-1]
       
    fac0 = (1.0*1000/press[0])**(1.0*rgas/cp)
    pottemp.append(1.0*abstemp[0]*fac0)
        
    height.append((1.0*rgas/ggr)*abstemp[0]*np.log(1.0*press0/press[0]))
                        
    for i in range(end-1):                
         height.append(height[i] + (1.0*rgas/ggr)*.5*(abstemp[i+1]+ abstemp[i])*np.log(1.0*press[i]/press[i+1]))
         fac = (1.0*1000/press[i+1])**(1.0*rgas/cp)
         pottemp.append(1.0*abstemp[i+1]*fac)             
    return np.array(pottemp), np.array(height)

def Get_Data(filename, start, stop):            
    """Takes output from txt file and converts to a np.array

    Argument:
    filename -- filename including end bit, enclosed in single inverted commas
    start, stop -- integers, lines to read if there's a header etc
    
    Returns:
    array -- np.array of file contents

    """ 
    file = open(filename)
    if start > 0:
        lines_to_read = file.readlines()[start:stop]
    else:
        lines_to_read = file    
    data = []    
    for line in lines_to_read:
        
        text = line.split()
        values = [float(i) for i in text]
        
        data.append(values)
    file.close()
    array = np.array(data)
    
    return array

def from_lmo():
     """Pulls output relating to Monin Obvukov Lenth
  
    Returns:
     -- array    
    
    """
     print 'Need to edit filepath for lmo.txt'
     txtfile_list = ["/tera/phil/nchaparr/sam_ensemble/sam_case" + str(i+2) + "/OUT_STAT/lmo.txt" for i in range(9)]
     array_list = []
     for txtfile in txtfile_list:
          array = np.genfromtxt(txtfile)
          [columns] = array.shape
          array = np.reshape(array, [1.0*columns/4, 4])
          array_list.append(array)
          
     array = Ensemble1_Average(array_list)     
     return array          

def Plot_nc(theFile, imax, theAx):
    """plots from nc file

    Arguments:
    theFile -- nc file
    imax -- integer, max index to loop through for plotting 
    
    """
    ncdata = Dataset(theFile,'r')
    press = 1.0*ncdata.variables['lev'][...]/100
    time = ncdata.variables['tsec'][...]
    abstemp = np.squeeze(ncdata.variables['T'][...]) #absolute temperature
    press0 = 1000    
           
    for i in range(len(time)):        
         if i<imax:            
            [pottemp, height] = abs2pot(press0, press, abstemp[i,:], 1) #to potential temperature          
            
            theAx.plot(pottemp, height,'o', label = 'from initial nc file')            
         


def Do_Plot(fignum, title, ylabel, xlabel, sub):
    """Returns subplot instance with title and axis lables

    Argument:
    fignum -- integer, identifies subplot instance
    title -- string, obvious
    ylable, xlabel -- obvious
    sub -- eg 111 for one single plot

    Returns:
    Ax -- subplot instance
    
    """
    Fig = plt.figure(fignum) 
    Fig.clf()
    Ax = Fig.add_subplot(sub) 
    Ax.set_title(title)
    Ax.set_xlabel(xlabel)
    Ax.set_ylabel(ylabel)
    return Ax

def Plot_Save(Ax, xvals, yvals, Yvals, ysymbol, Ysymbol, firsttag, secondtag, figname):
    """Plots two sets of y values against one set of x values and saves the figure. 

    Argument:
    Ax -- subplot instance name    
    xvals, yvals, Yvals -- obvious
    ysymbol, Ysymbol -- matplotlib symbol enclosed in single inverted commas (i think)
    firsttag, secondtag -- legends
    figname -- name to be saved as   
    
    """
    Ax.plot(xvals, yvals, ysymbol, label = firsttag)
    Ax.plot(xvals, Yvals, Ysymbol, label = secondtag)
    plt.legend(loc = 'upper left', prop={'size':8})    
    plt.show()
    plt.savefig(figname)    


def indirect_sort(list1, list2):
    """

    Argument:
    filename -- filename including end bit, enclosed in single inverted commas

    Returns:
    array -- np.array of file contents

    """
    sort_index=np.argsort(list1)
    new_list2 = []
    for i in range(len(list1)):
        
        new_list2.append(list2[sort_index[i]])
    return new_list2

def Horizontal_Average(array):
     """Gets horizontal average of 64, , array

    Arguments:
    array -- 64,  array

    Returns:
    array_bar -- 64 array    
    
    """
     array_bar = []
     
     [zrows, ycols, xcols] = array.shape
     for i in range(zrows):
          avvals = array[i, :, :]          
          avvals = np.mean(avvals, axis=1)     
          avvals = np.mean(avvals, axis=0)     
          array_bar.append(avvals)
         
     return np.array(array_bar)

def Domain_Grad(array, height):
     """Gets vertical gradients and height of max gradient 

    Arguments:
    --  array

    Returns:
     -- array    
    
    """
     from scipy import signal
     height3d = np.zeros_like(array)
     [rows, ycols, xcols] = array.shape
     for i in range(rows):
          height3d[i,:,:] = height[i]  
          
     
     dvar = np.gradient(array)[0]
     
     dheight = np.gradient(height3d)[0]
     dvardz = np.divide(dvar, dheight)
        
        
     grad_peaks = np.zeros((ycols, xcols))
     bot_index = np.where(abs(100 - height) < 26.)[0][0]
     
     steps = np.zeros_like(array)
     count = 0
     for i in range(ycols):
          
          for j in range(xcols):
               
               max_grad = np.amax(dvardz[bot_index:, i, j])
               
               index = np.where(dvardz[:, i, j] - max_grad == 0)[0][0]
               
               BLheight = height[index]               
               grad_peaks[i, j] = BLheight
               count = count+1
               
     return dvardz, grad_peaks

def Bin_Peaks(peaks, heights):
     """Puts gradient peak heights into bins 

    Arguments:
    --  

    Returns:
     --     
    
    """
     bin_vols=np.zeros_like(heights)
     [xpos, ypos] = peaks.shape
     print xpos, ypos
     for i in range(xpos):
          for j in range(ypos):
                            
               index = np.where(heights - peaks[i, j] == 0)[0][0]
               
               bin_vols[index] = bin_vols[index] + 1               
     return bin_vols          


def Ensemble_Average(list):
     """Gets enseble average of a list of arrays

    Arguments:
    list -- list of 64, 64, 64 arrays

    Returns:
    ens_avs -- 64, 64, 64 array
    
    """     
     ens_avs = np.zeros([128, 64, 64])
     for i in range(128): #vertical
          for j in range(64):#horizontal 
               for k in range(64): #horizontal y
                    vals = []                              
                    for l in range(len(list)):
                         val = list[l][i, j, k]#i, j, kth elemement from var array l                                        
                         vals.append(val)                
                    avval = 1.0*sum(vals)/len(vals) #average over l arrays                             
                    ens_avs[i, j, k] = avval
     return ens_avs

def Ensemble1_Average(list):
    """Gets enseble average of a list of arrays

    Arguments:
    list -- list of arrays

    Returns:
    ens_avs -- array

    """
    to_av = list[0]
        
    for k in range(len(list)-1):
         #print k, 'array sizes', to_av.shape, list[k+1].shape 
         to_av = np.add(to_av, list[k+1])
    ens_avs = 1.0*to_av/len(list)
    
    return ens_avs


def pert_h_Average(list):
    """Gets average of pos or neg perts at each height level from a list of arrays

    Arguments:
    list -- list of arrays

    Returns:
    avs -- array

    """
    A = list[0]
    print A.shape, "shape of array to be averaged"
    A = np.ma.average(A, axis=1, weights=A.astype(bool))
    print A.shape, "shape of array to be averaged"
    to_av = np.ma.average(A, axis=1, weights=A.astype(bool))
    print to_av.shape, "shape of array to be averaged"
    for k in range(len(list)-1):
         #print k, 'array sizes', to_av.shape, list[k+1].shape
         B = list[k+1]
         B = np.ma.average(B, axis=1, weights=B.astype(bool))
         B = np.ma.average(B, axis=1, weights=B.astype(bool))
         to_av = np.add(to_av, B)
    avs = 1.0*to_av/len(list)
    print avs.shape
        
    return avs


def Flux_Quad_Slow(wpert, thetapert):
    """
    Separates fluxes into quadratns
    
    Arguments:
    wpert -- array of w perturbations
    thetapert -- array of theta perturbations

    Returns:
    up_warm, down_warm, up_cold, down_cold -- arrays, np.nans are fillers  

    """
    
    [rows, columnsx, columnsy] = wpert.shape
    array_list = [np.zeros_like(wpert) for i in range(4)]
    for j in range(rows):
         for k in range(columsx):
              for l in range(columsy):
                   if wpert[j,k,l]>0 and thetapert[j,k,l]>0:
                        array_list[0][i,j,k], array_list[1][i,j,k], array_list[2][i,j,k], array_list[3][i,j,k] = wpert[j,k,l]*thetapert[j,k,l], np.nan, np.nan, np.nan
                   elif wpert[j,k,l]<0 and thetapert[j,k,l]>0:
                        array_list[0][i,j,k], array_list[1][i,j,k], array_list[2][i,j,k], array_list[3][i,j,k] = np.nan, wpert[j,k,l]*thetapert[j,k,l], np.nan, np.nan
                   elif wpert[j,k,l]>0 and thetapert[j,k,l]<0:
                        array_list[0][i,j,k], array_list[1][i,j,k], array_list[2][i,j,k], array_list[3][i,j,k] = np.nan, np.nan, wpert[j,k,l]*thetapert[j,k,l], np.nan
                   else:
                       array_list[0][i,j,k], array_list[1][i,j,k], array_list[2][i,j,k], array_list[3][i,j,k] = np.nan, np.nan, np.nan, wpert[j,k,l]*thetapert[j,k,l]                    

    [avup_warm, avdown_warm, avup_cold, avdown_cold] = [Horizontal_Average(array) for array in array_list ]                    
    
                                               
    return [avup_warm, avdown_warm, avup_cold, avdown_cold]



def Get_Var_Arrays(ncfolder, ncfilename, dump_time, case_number):
    #TODO: make more modular, eg add option to dump txt files
    #TODO: combine functions with class Get_Var_Arrays1
     """Pulls output from an ensemble cases

    Arguments:
    dump_time, case_number -- time of output eg '0000000720', obvious eg 1

    Returns:
    var_bar -- 64 array of horizontally averaged, ensemble averages or perturbations (covariances)
    
    """
     #create list of filenames for given dump_time
     ncfile = ncfolder + str(case_number) + ncfilename + dump_time + ".nc"
     
     #create lists for variable arrays from each case
     
     thefile = ncfile
     print thefile     
     ncdata = Dataset(thefile,'r')
     wvel = np.squeeze(ncdata.variables['W'][...])
     
     press = np.squeeze(ncdata.variables['p'][...])#pressure already horizontally averaged
     height = np.squeeze(ncdata.variables['z'][...])
     temp = np.squeeze(ncdata.variables['TABS'][...])
     #tracer = np.squeeze(ncdata.variables['TRACER'][...])
     ncdata.close()
     
          #calculate thetas
     theta = np.zeros_like(temp)
     thetafact = np.array([(1.0*1000/k)**(1.0*287/1004) for k in press])
     [zvals, yvals, xvals] = theta.shape
     for j in range(zvals): #TODO: do a theta.shape, to the z dimension
         theta[j, :, :] = temp[j, :, :]*thetafact[j]     
     
     return wvel, theta, theta, height


def Flux_Quad(wpert, thetapert):
    """
    Separates fluxes into quadrants
    
    Arguments:
    wpert -- array of w perturbations
    thetapert -- array of theta perturbations

    Returns:
    [up_warm, down_warm, up_cold, down_cold] -- arrays, np.nans are fillers  

    """
    
    [rows, columnsx, columnsy] = wpert.shape
    [up, down, warm, cold] = [np.zeros_like(wpert) for i in range(4)]
    up[wpert>0], down[wpert<0], warm[thetapert>0], cold[thetapert<0] = wpert[wpert>0], wpert[wpert<0], thetapert[thetapert>0], thetapert[thetapert<0]
    upwarm, downwarm, upcold, downcold = np.multiply(up, warm), np.multiply(down, warm), np.multiply(up, cold), np.multiply(down, cold)
                                                             
    return [upwarm, downwarm, upcold, downcold]

def Flux_Quad_Thetas(wpert, thetapert):
    """
    Separates fluxes into quadrants
    
    Arguments:
    wpert -- array of w perturbations
    thetapert -- array of theta perturbations

    Returns:
    [up_warm, down_warm, up_cold, down_cold] -- arrays, np.nans are fillers  

    """
    
    [rows, columnsx, columnsy] = wpert.shape
    [up, down, warm, cold] = [np.zeros_like(wpert) for i in range(4)]
    up[wpert>0], down[wpert<0], warm[thetapert>0], cold[thetapert<0] = wpert[wpert>0], wpert[wpert<0], thetapert[thetapert>0], thetapert[thetapert<0]
    up[up>0], down[down<0] = 1, 1    
    upwarm, downwarm, upcold, downcold = np.multiply(up, warm), np.multiply(down, warm), np.multiply(up, cold), np.multiply(down, cold)
                                                             
    return [upwarm, downwarm, upcold, downcold]

def Flux_Quad_Wvels(wpert, thetapert):
    """
    Separates fluxes into quadrants
    
    Arguments:
    wpert -- array of w perturbations
    thetapert -- array of theta perturbations

    Returns:
    [up_warm, down_warm, up_cold, down_cold] -- arrays, np.nans are fillers  

    """
    
    [rows, columnsx, columnsy] = wpert.shape
    [up, down, warm, cold] = [np.zeros_like(wpert) for i in range(4)]
    up[wpert>0], down[wpert<0], warm[thetapert>0], cold[thetapert<0] = wpert[wpert>0], wpert[wpert<0], thetapert[thetapert>0], thetapert[thetapert<0]
    warm[warm>0], cold[cold<0] = 1, 1    
    upwarm, downwarm, upcold, downcold = np.multiply(up, warm), np.multiply(down, warm), np.multiply(up, cold), np.multiply(down, cold)
                                                             
    return [upwarm, downwarm, upcold, downcold]


def Get_CBLHeights(heights, press, thetas, theta0, wvelthetapert, gamma, flux_s, top_index):
    """
    Gets heights based on dthetdz and flux
    
    Arguments:
    wpert -- array of w perturbations
    thetapert -- array of theta perturbations

    Returns:
    [up_warm, down_warm, up_cold, down_cold] -- arrays, np.nans are fillers  

    """
    
    dheight = np.diff(heights)
    dtheta = np.diff(thetas)
    dthetadz=np.divide(dtheta, dheight)
    element0=np.array([0])
    dthetadz=np.hstack((element0, dthetadz))*1.0/gamma
    rhow=calc_rhow(press, heights, thetas[0])
    fluxes=np.multiply(wvelthetapert, rhow)*1004.0/flux_s
    
    #where gradient is greater than zero    
    for j in range(len(dthetadz[:top_index])-1):
        if (dthetadz[j+1] >.03) and (dthetadz[j] >= 0):
            dtheta_index_b = j+1
            
            break

    #where gradient resumes as gamma
    dtheta_index_t = 999
    for k in range(len(dthetadz[:top_index])-1):
        print k
        print ""
        print dthetadz[k-1], dthetadz[k+1], dthetadz[k+2]
        print ""
        print np.abs(dthetadz[k+1]-1), np.abs(dthetadz[k+2]-1)
        if np.abs(dthetadz[k+2]-1)<.03 and np.abs(dthetadz[k+1]-1)<.03 and dthetadz[k-1]>1:            
            dtheta_index_t = k+1
            
            break

    #Hacky fix for when the upper theta gradient profiles are wonky
    if dtheta_index_t == 999:
         for k in range(len(dthetadz[:top_index])-1):
           #   print dthetadz[k-1], dthetadz[k+1], dthetadz[k+2]
           #   print ""
           #   print np.abs(dthetadz[k+1]-1), np.abs(dthetadz[k+2]-1)
              if np.abs(dthetadz[k+2]-1)<.04 and np.abs(dthetadz[k+1]-1)<.04 and dthetadz[k-1]>1:            
                   dtheta_index_t = k+1                        
                   break
         
         
    #now fluxes    
    
    for l in range(len(dthetadz)-1):
        if (fluxes[l+1] <= .0) and (fluxes[l] > 0):
            flux_index_b = l+1
            break
        
    for m in range(len(dthetadz[0:top_index])-1):
         #print fluxes[m+1], fluxes[m], fluxes[m-1]
         if (abs(fluxes[m+1]) < 0.01) and (abs(fluxes[m+2]) < 0.01) and (fluxes[m] < 0) and (fluxes[m-1] < 0):
            flux_index_t = m+1
            break
    #print flux_index_t
    
    eltop_dthetadz = heights[dtheta_index_t]
    elbot_dthetadz = heights[dtheta_index_b]

    eltop_flux = heights[flux_index_t]
    elbot_flux = heights[flux_index_b]    

    h = heights[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]]
    h_flux = heights[np.where(wvelthetapert - np.amin(wvelthetapert) == 0)[0][0]]
    
    deltatheta_f = theta0[np.where(wvelthetapert - np.amin(wvelthetapert) == 0)[0][0]] - thetas[np.where(wvelthetapert - np.amin(wvelthetapert) == 0)[0][0]]
    deltatheta = theta0[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]] - thetas[np.where(dthetadz[0:top_index] - np.amax(dthetadz[0:top_index]) == 0)[0][0]]
    Deltatheta_f = thetas[flux_index_t] - thetas[flux_index_b]
    Deltatheta = thetas[dtheta_index_t] - thetas[dtheta_index_b]

    mltheta = np.mean(thetas[0:dtheta_index_b])
    
    return [elbot_dthetadz, h, eltop_dthetadz, elbot_flux ,h_flux  ,eltop_flux, Deltatheta, Deltatheta_f, deltatheta, deltatheta_f, mltheta]

    
    
    
