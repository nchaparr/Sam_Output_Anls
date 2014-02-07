from netCDF4 import Dataset
import numpy as np

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
          return wvelperts_list

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

    
    
    
