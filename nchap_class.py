from netCDF4 import Dataset
import numpy as np
import nchap_fun as nc
            

class For_Plots:
     """
        for pulling data for plot_height.py ie plots representing all runs on one ax instance
     
     """
     def __init__(self, Run_Date):
          self.Run_Date = Run_Date          
          self.path = "/newtera/tera/phil/nchaparr/python/Plotting/" + Run_Date + "/data/"

     def get_file(self, dump_time, filename):
          the_file = self.path + filename + dump_time
          return the_file

     def save_file(self, array, filename):
          np.savetxt(self.path + filename, array, delimiter=' ')
          
     def dhdtplot(self):
          dhdtplot = np.genfromtxt(self.path + "dhdtinvriplt.txt")
          return dhdtplot

     def HistVars(self):
          HistVars = np.genfromtxt(self.path + "ml_height_hist_vars")
          return HistVars

     def AvProfVars(self):
          AvProfVars = np.genfromtxt(self.path + "AvProfLims")
          return AvProfVars

     def rinovals(self):
          rinovals = np.genfromtxt(self.path + "invrinos")
          return rinovals

     def scaled_we_plot(self):
          scaled_we_plot = np.genfromtxt(self.path + "dhdtinvriplt.txt")
          return scaled_we_plot
 

     def Deltah(self):
          AvProfVars = np.genfromtxt(self.path + "AvProfLims")
          Deltah = np.subtract(AvProfVars[:,2], AvProfVars[:,0])
          #Deltah0 = np.divide(Deltah0, AvProfVars[:,1])
          return Deltah
          
     def Deltah_over_h(self):
          AvProfVars = np.genfromtxt(self.path + "AvProfLims")
          Deltah = np.subtract(AvProfVars[:,2], AvProfVars[:,0])
          Deltah_over_h = np.divide(Deltah, AvProfVars[:,1])
          return Deltah_over_h         

     def Get_and_save_dhdt(self, Times, Heights, Wstars, Invrinos):
         import numpy as np
         FitFunc=np.polyfit(Times, Heights, 2, full=False)
         Fit = FitFunc[0]*Times**2 + FitFunc[1]*Times + FitFunc[2]
         dhdt =1.0*(2*FitFunc[0]*Times + FitFunc[1])/3600
         scaled_dhdt = np.divide(dhdt, Wstars)
         dhdtinvriplt = np.vstack((Invrinos, scaled_dhdt))
         np.savetxt(self.path+"dhdtinvriplt.txt", dhdtinvriplt, delimiter=' ')


class Get_Var_Arrays1:
     """
        for pulling velociy perturpations, temperature and pressure
        from nc files from the ensemble
       
     """
     def __init__(self, path1, path2, dump_time):
          self.path1 = path1
          self.path2 = path2          
          self.dump_time = dump_time
          self.nc_file_list = [path1 + str(i+1) + path2 + dump_time + ".nc" for i in range(10)]
     
     
     def get_wvelperts(self):
          wvelperts_list = []
          for i in range(len(self.nc_file_list)): #loop over list of nc files, not efficient to do this for each variable but can't think of better way right now
               thefile = self.nc_file_list[i]
               print thefile
               ncdata = Dataset(thefile,'r')          
               wvelperts_list.append(np.squeeze(ncdata.variables['W'][...]))
               ncdata.close()
          return wvelperts_list
     

     def get_uvelperts(self):
          uvelperts_list = []
          for i in range(len(self.nc_file_list)): #loop over list of nc files, not efficient to do this for each variable but can't think of better way right now
               thefile = self.nc_file_list[i]
               print thefile
               ncdata = Dataset(thefile,'r')          
               uvelperts_list.append(np.squeeze(ncdata.variables['U'][...]))
               ncdata.close()
          return uvelperts_list

     def get_vvelperts(self):
          vvelperts_list = []
          for i in range(len(self.nc_file_list)): #loop over list of nc files, not efficient to do this for each variable but can't think of better way right now
               thefile = self.nc_file_list[i]
               print thefile
               ncdata = Dataset(thefile,'r')          
               vvelperts_list.append(np.squeeze(ncdata.variables['V'][...]))
               ncdata.close()
          return vvelperts_list


     def get_thetas(self):
          thetas_list = []
          press_list = []
          for i in range(len(self.nc_file_list)): #loop over list of nc files, not efficient to do this for each variable but can't think of better way right now
               thefile = self.nc_file_list[i]
               print thefile
               ncdata = Dataset(thefile,'r')          
               temp = np.squeeze(ncdata.variables['TABS'][...])
               press = np.squeeze(ncdata.variables['p'][...])
               [znum, ynum, xnum] = temp.shape
               theta = np.zeros_like(temp) #TODO: make this into a function
               thetafact = np.array([(1.0*1000/k)**(1.0*287/1004) for k in press])
               for j in range(znum):
                    theta[j, :, :] = temp[j, :, :]*thetafact[j]
               thetas_list.append(theta)
               press_list.append(press)
               ncdata.close()
          return thetas_list, press_list

     def get_height(self):
          thefile = self.nc_file_list[0]
          print thefile
          ncdata = Dataset(thefile,'r')          
          height = np.squeeze(ncdata.variables['z'][...])
          ncdata.close()
          return height

     def get_wvelthetaperts(self):
          wvels_list = self.get_wvelperts()
          thetas_list, press_list = self.get_thetas()
          #print 'checking thetas_list', len(thetas_list), thetas_list[0].shape
          ens_avthetas = nc.Ensemble1_Average(thetas_list)
          wvelthetaperts_list = []
          for i in range(len(wvels_list)):  
               [znum, ynum, xnum] = wvels_list[i].shape
               thetapert_rough = np.subtract(thetas_list[i], ens_avthetas)
               thetapert = np.zeros_like(thetapert_rough)
               for j in range(znum):#something like this is done in statistics.f90, staggered grid!
                    if j == 0:
                         thetapert[j,:,:] = thetapert_rough[j,:,:]
                    else:
                         thetapert[j,:,:] = 0.5*np.add(thetapert_rough[j,:,:], thetapert_rough[j-1,:,:])
               wvelpert = wvels_list[i]     
               wvelthetapert = np.multiply(wvelpert, thetapert)
               wvelthetaperts_list.append(wvelthetapert)

          return wvelthetaperts_list

     def get_thetaperts(self):
          #wvels_list = self.get_wvelperts()
          thetas_list, press_list = self.get_thetas()
          #print 'checking thetas_list', len(thetas_list), thetas_list[0].shape
          ens_avthetas = nc.Ensemble1_Average(thetas_list)
          thetaperts_list = []
          for i in range(len(thetas_list)):  
               [znum, ynum, xnum] = thetas_list[i].shape
               thetapert_rough = np.subtract(thetas_list[i], ens_avthetas)
               thetapert = np.zeros_like(thetapert_rough)
               for j in range(znum):#something like this is done in statistics.f90, staggered grid!
                    if j == 0:
                         thetapert[j,:,:] = thetapert_rough[j,:,:]
                    else:
                         thetapert[j,:,:] = 0.5*np.add(thetapert_rough[j,:,:], thetapert_rough[j-1,:,:])
               #wvelpert = wvels_list[i]     
               #wvelthetapert = np.multiply(wvelpert, thetapert)
               thetaperts_list.append(thetapert)

          return thetaperts_list
     
     def get_sqvel(self, vel_dir):
          if vel_dir == 'w':
               vel_list = self.get_wvelperts()
          elif vel_dir == 'v':
               vel_list = self.get_vvelperts()
          elif vel_dir == 'u':
               vel_list = self.get_uvelperts()
               
          velpertsq_list = []

          for i in range(len(vel_list)):
               velpert = vel_list[i]
               velpertsq = np.multiply(velpert, velpert)
               velpertsq_list.append(velpertsq)

          return velpertsq_list
    
