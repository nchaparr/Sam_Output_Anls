import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

run_key={'Nov302013':{'params':(100,5)},         
         'Dec142013':{'params':(100,10)},        
         'Dec202013':{'params':(60,5)},          
         'Dec252013':{'params':(60,2.5)},        
         'Jan152014_1':{'params':(150,5)},         
         'Mar12014':{'params':(60,10)},         
         'Mar52014':{'params':(150,10)}}        

def gm_vars(h,surface_flux,gamma):
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
    Re0=(L0*B0)**(1./3.)
    wstar=(g*h/theta_0*flux)**(1./3)
    thetastar=flux/wstar
    wstar_gm=(B0*h)**(1./3.)
    c_gamma=0.55
    delta=c_gamma*wstar/N
    return L0,N,B0,delta

def find_height(height_m,time_sec,N,L0):
    zenc=L0*(2*time_sec*N)**0.5
    height_nd=height_m/zenc
    return height_nd,zenc

if __name__ == "__main__":
    case_list=[]
    datadir='/tera/phil/nchaparr/python/Plotting'
    h=500.
    for case,run_dict in run_key.items():
        if case=='Nov302013':
            time_int=900
        else:
            time_int=600
        surface_flux,gamma=run_dict['params']
        L0,N,B0,delta=gm_vars(h,surface_flux,gamma)
        filename='{}/{}/data/AvProfLims'.format(datadir,case)
        out=np.genfromtxt(filename)
        columns=['h0','h','h1','zf0','zf','zf1','deltatheta','mltheta']
        df=pd.DataFrame.from_records(out,columns=columns)
        time_end=28800 
        num_times=len(df)
        time_beg=time_end -time_int*num_times
        time_sec=np.linspace(time_beg,time_end,num_times)
        run_dict['df']=df
        run_dict['L0']=L0
        run_dict['N']=N
        run_dict['df']['time_sec']=time_sec
        run_dict['df']['time_nd']=time_sec*N
        case_list.append((case,L0))
    case_list.sort(key=lambda case: case[1])
    first=case_list[0][0]
    L0,N=run_key[first]['L0'],run_key[first]['N']
    df=run_key[first]['df']
    height_nd,zenc=find_height(df['h'],df['time_sec'],N,L0)
    h0_nd,zenc=find_height(df['h0'],df['time_sec'],N,L0)
    h1_nd,zenc=find_height(df['h1'],df['time_sec'],N,L0)
    df['delhtop']=(df['h1'] - df['h'])/df['h']
    df['delhbot']=(df['h'] - df['h0'])/df['h']
    df['delhtot']=(df['h1'] - df['h0'])/df['h']
    df['delgm']=0.55*(zenc/L0)**(-2/3.)  #eq. 26
    df['zenc']=zenc
    plt.close('all')
    fig_h,ax_h=plt.subplots(1,1)
    fig_h0,ax_h0=plt.subplots(1,1)
    fig_h1,ax_h1=plt.subplots(1,1)
    fig_delhtop,ax_delhtop=plt.subplots(1,1)
    fig_delgm,ax_delgm=plt.subplots(1,1)
    fig_delhbot,ax_delhbot=plt.subplots(1,1)
    fig_delhtot,ax_delhtot=plt.subplots(1,1)
    fig_delhtot_rt,ax_delhtot_rt=plt.subplots(1,1)
    label='{:3.1f}'.format(L0)
    ax_h.plot(df['time_nd'],height_nd,label=label)
    ax_h0.plot(df['time_nd'],h0_nd,label=label)
    ax_h1.plot(df['time_nd'],h1_nd,label=label)
    ax_delhtop.plot(df['time_nd'],df['delhtop'],label=label)
    ax_delhbot.plot(df['time_nd'],df['delhbot'],label=label)
    ax_delhtot.plot(df['time_nd'],df['delhtot'],label=label)
    ax_delgm.plot(df['time_nd'],df['delgm'],label=label)
    ax_delhtot_rt.plot(df['time_sec'],df['delhtot'],label=label)
    for casename,L0 in case_list[1:]:
        N=run_key[casename]['N']
        df=run_key[casename]['df']
        height_nd,zenc=find_height(df['h'],df['time_sec'],N,L0)
        df['zenc']=zenc
        h0_nd,zenc=find_height(df['h0'],df['time_sec'],N,L0)
        h1_nd,zenc=find_height(df['h1'],df['time_sec'],N,L0)
        df['delhtop']=(df['h1'] - df['h'])/df['h']
        df['delhbot']=(df['h'] - df['h0'])/df['h']
        df['delhtot']=(df['h1'] - df['h0'])/df['h']
        df['delgm']=0.55*(zenc/L0)**(-2/3.)  #eq. 26
        label='{:3.1f}'.format(L0)
        ax_h.plot(df['time_nd'],height_nd,label=label)
        ax_h0.plot(df['time_nd'],h0_nd,label=label)
        ax_h1.plot(df['time_nd'],h1_nd,label=label)
        ax_delhtop.plot(df['time_nd'],df['delhtop'],label=label)
        ax_delhbot.plot(df['time_nd'],df['delhbot'],label=label)
        ax_delhtot.plot(df['time_nd'],df['delhtot'],label=label)
        ax_delhtot_rt.plot(df['time_sec'],df['delhtot'],label=label)

    plot_dict=dict(h=ax_h,h1=ax_h1,h0=ax_h0,
                   delhtop=ax_delhtop,delhbot=ax_delhbot,delhtot=ax_delhtot,delgm=ax_delgm,
                   delhtot_rt=ax_delhtot_rt)
    titles=dict(h='non-dimensional h',h1='non-dimensional h1',h0='non-dimensional h0',
                delhtop='(h1 - h)/h',delhbot='(h - h0)/h',delhtot='(h1 - h0)/h',
                delgm='$\delta/z_{enc}$',delhtot_rt='(h1 - h0)/h vs. dimensional time')
    ylims=dict(h=(0.6,1.5),h1=(0.6,1.5),h0=(0.6,1.5),
               delgm=(0.0,0.7),delhtop=(0.0,0.7),delhbot=(0.0,0.7),
               delhtot=(0.0,0.7),delhtot_rt=(0.0,0.7))
    for key,the_ax in plot_dict.items():
        the_ax.set_ylim(ylims[key])
        the_ax.set_title(titles[key])
        the_ax.legend(loc='best')
        if key == 'delhtot_rt':
            the_ax.set_xlabel('time (sec)')
        else:
            the_ax.set_xlabel('time*N')
    cases=[name[0] for name in case_list]
    df_cases=pd.DataFrame(cases,columns=['name'])
    #
    #  params keyword is (flux,gamma) tuple
    #
    gammas=[run_key[case]['params'][1] for case in cases]
    fluxes=[run_key[case]['params'][0] for case in cases]
    l0=[run_key[case]['L0'] for case in cases]
    N=[run_key[case]['N'] for case in cases]
    period=[1./Nval/60. for Nval in N]   #minutes
    df_cases['L0']=l0
    df_cases['Nperiod']=period
    df_cases['fluxes']=fluxes
    df_cases['gammas']=gammas
    for key,axis in plot_dict.items():
        filename='{}.png'.format(key)
        axis.figure.savefig(filename)
    plt.show()
    with pd.HDFStore('paper_table.h5','w') as store:
        store.put('cases',df_cases,format='table')

        

        



#Files = /newtera/tera/phil/nchaparr/python/Plotting/rundate/data/AvProfLims
#indices of h0, h, h1 = 0, 1, 2
#rundate = [Dec142013, Nov302013, Dec202013, Dec252013, Jan152014_1, Mar12014, Mar52014]
#sfcflx/gamma = [100/10, 100/5, 60/5, 60/2.5, 150/5, 60/10, 150/10]
#time increments for all except Nov302013 = 600, for Nov302013 = 900
#end time = 28800

