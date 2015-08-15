import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from collections import OrderedDict as od
h5new='good.h5'
#
# get the root attributes
#

df_dict=od()
with pd.HDFStore(h5new,'r') as store:
    df_overview=store.get('/df_overview')
    varnames=list(store.get('/var_names'))
    case_names=list(store.get('/date_names'))
    time600=np.array(store.get('/time600'))
    time900=np.array(store.get('/time900'))
    for case in case_names:
        for name in ['AvProfLims']:
            nodename='{}/{}'.format(case,name)
            df_dict[case,name]=store.get(nodename)
#
# get N and L0 for each case
#
Nvals=df_overview['N']
names=df_overview['name']
L0=df_overview['L0']

run_key={}
for item in df_overview.to_dict('records'):
    run_key[item['name']]=dict(params=(item['fluxes'],item['gammas']))

def gm_vars(surface_flux,gamma):
    rho=1.
    cp=1004.
    g=9.8
    flux=surface_flux/(rho*cp)  #from W/m^2 to m K/s
    theta_0=300.  #K
    B0=flux*g/theta_0
    gamma=gamma/1000.  #K/m
    g=9.8  #m/s^2
    N2=g/theta_0*gamma  #s**(-2)
    N=N2**0.5
    L0=(B0/N**3.)**0.5  #gm eqn 3
    #Re0=(L0*B0)**(1./3.)
    # wstar=(g*h/theta_0*flux)**(1./3)
    # c_gamma=0.55
    # delta=c_gamma*wstar/N
    #thetastar=flux/wstar
    #wstar_gm=(B0*h)**(1./3.)
    return L0,N,B0

if __name__ == "__main__":
    case_list=[]
    for case in run_key.keys():
        run_dict=run_key[case]
        surface_flux,gamma=run_dict['params']
        L0,N,B0=gm_vars(surface_flux,gamma)
        df=df_dict[case,'AvProfLims']
        run_dict['df']=df
        run_dict['L0']=L0
        run_dict['N']=N
        run_dict['df']['time_secs']=df['times']
        run_dict['df']['time_nd']=df['Ntimes']
        run_dict['df']['h0_nd']=df['h0']/df['zenc']
        case_list.append((case,L0))

    case_list.sort(key=lambda case: case[1])
    plt.close('all')
    plotlist=['h','h1','h0','delhtop','delhbot','delhtot','delgm',
              'delhtot_rt','zf','zf0','zf1','delzfbot']
    xy_dict={'h':('time_nd','h_nd'),'h1':('time_nd','h1_nd'),'h0':('time_nd','h0_nd'),'delhtop':('time_nd','delhtop'),
             'delhbot':('time_nd','delhbot'),'delhtot':('time_nd','delhtot'),'delgm':('time_nd','delgm'),
             'delhtot_rt':('time_secs','delhtot_rt'),'zf':('time_nd','zf_nd'),
             'zf0':('time_nd','zf0_nd'),'zf1':('time_nd','zf1_nd'),'delzfbot':('time_nd','delzfbot')}
    titles=dict(h='non-dimensional h',h1='non-dimensional h1',h0='non-dimensional h0',
                delhtop='(h1 - h)/h',delhbot='(h - h0)/h',delhtot='(h1 - h0)/h',
                delgm='$\delta/z_{enc}$',delhtot_rt='(h1 - h0)/h vs. dimensional time',
                zf='non-dimensional zf',zf0='non-dimensional zf0',zf1='non-dimensional zf1',
                delzfbot='(zf - zf0)/zf')
    ylims=dict(h=(0.6,1.5),h1=(0.6,1.5),h0=(0.6,1.5),zf=(0.6,1.5),zf1=(0.6,1.5),zf0=(0.6,1.5),
               delgm=(0.0,0.7),delhtop=(0.0,0.7),delhbot=(0.0,0.7),
               delhtot=(0.0,0.7),delhtot_rt=(0.0,0.7),delzfbot=(0,0.7))
    plot_dict={}
    for plot in plotlist:
        fig,ax=plt.subplots(1,1)
        plot_dict[plot]=ax

    for casename,L0 in case_list:
        L0,N=run_key[casename]['L0'],run_key[casename]['N']
        df=df_dict[casename,'AvProfLims']
        for key in ['h','h0','h1','zf','zf0','zf1']:
            nd_key='{}_nd'.format(key)
            df[nd_key]=df[key]/df['zenc']
        df['delhtop']=(df['h1'] - df['h'])/df['h']
        df['delhbot']=(df['h'] - df['h0'])/df['h']
        df['delhtot']=(df['h1'] - df['h0'])/df['h']
        df['delgm']=0.55*(df['zenc']/L0)**(-2/3.)  #eq. 26
        df['delhtot_rt']=(df['h1'] - df['h0'])/df['h']
        df['delzfbot']=(df['zf'] - df['zf0'])/df['zf']
        label='{:3.1f}'.format(L0)
        for plot in plotlist:
            xvals=df[xy_dict[plot][0]]
            yvals=df[xy_dict[plot][1]]
            plot_dict[plot].plot(xvals,yvals,label=label)
    for plot,the_ax in plot_dict.items():
        the_ax.set_ylim(ylims[plot])
        the_ax.set_title(titles[plot])
        the_ax.legend(loc='best')
        if plot == 'delhtot_rt':
            the_ax.set_xlabel('time (sec)')
        else:
            the_ax.set_xlabel('time*N')
    for key,axis in plot_dict.items():
        filename='{}.png'.format(key)
        axis.figure.savefig(filename)
    plt.show()
        

