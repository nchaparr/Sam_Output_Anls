import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import cm
from matplotlib.colors import Normalize
#from panormalize import clippedNorm
import matplotlib
import numpy as np

import matplotlib.pyplot as plt
theFig=plt.figure(1)
theFig.clf()
theAx=theFig.add_subplot(121)
the_file='/tera/phil/nchaparr/sam_ensemble/sam_case2/OUT_STAT/NCHAPP1_testing_doscamiopdata.nc'
the_nc=Dataset(the_file)
tke=the_nc.variables['THETAV'][...]
#tke=np.ma.masked_array(tke,mask=the_mask)
z=the_nc.variables['z'][...]
np.savetxt('/tera/phil/nchaparr/python/Plotting/Sep302013/data/heights', z, delimiter=' ')
print z
time=the_nc.variables['time'][...]
nrows,cols=tke.shape
print nrows
print cols
for nrow in range(nrows):
    print nrow
    if np.mod(nrow, 4) == 0:
        #if nrow < 18:
        theAx.plot(tke[nrow,:],z)
theAx.set_title('THETAV')
plt.ylim(0, 2500)
plt.xlim(295, 305)


theAx=theFig.add_subplot(122)
the_file='/tera/phil/nchaparr/sam_ensemble/sam_case2/OUT_STAT/NCHAPP1_testing_doscamiopdata.nc'

the_nc=Dataset(the_file)
tke=the_nc.variables['TVFLUX'][...]
#tke=np.ma.masked_array(tke,mask=the_mask)
z=the_nc.variables['z'][...]
time=the_nc.variables['time'][...]
#print time
nrows,cols=tke.shape
for nrow in range(nrows):
    if np.mod(nrow, 4) == 0:
        #if nrow < 18:
            print nrow
            theAx.plot(tke[nrow,:],z)
theAx.set_title('TVFLUX')
plt.ylim(0, 2500)
theFig.canvas.draw()
#theFig.savefig('test8.png',dpi=150)

plt.show()
 

