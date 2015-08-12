import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

with pd.HDFStore('mar12014.h5','r') as store:
    nodename='/Mar12014/flux_prof'
    df_flux=store.get(nodename)
    nodename='/Mar12014/h0'
    df_h0=store.get(nodename)

plt.close('all')
hit=np.where(df_h0['times']==27600)
h0=float(df_h0.loc[hit]['h0'])
flux=df_flux.as_matrix([27600]).squeeze()
heights=df_flux.as_matrix(['height'])/h0
fig,ax=plt.subplots(1,1)
ax.plot(flux,heights)
ax.set_ylim(0,1.3)
plt.show()
