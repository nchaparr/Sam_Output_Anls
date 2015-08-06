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

def find_zenc(time_sec,N,L0):
    zenc=L0*(2*time_sec*N)**0.5
    return zenc

if __name__ == "__main__":
    case_list=[]
    datadir='/tera/phil/nchaparr/python/Plotting'
    columns=['h0','h','h1','zf0','zf','zf1','deltatheta','mltheta']
    for case,run_dict in run_key.items():
        if case=='Nov302013':
            time_int=900
        else:
            time_int=600
        surface_flux,gamma=run_dict['params']
        L0,N,B0=gm_vars(surface_flux,gamma)
        filename='{}/{}/data/AvProfLims'.format(datadir,case)
        out=np.genfromtxt(filename)
        df=pd.DataFrame.from_records(out,columns=columns)
        time_end=28800 
        num_times=len(df)
        time_beg=time_end -time_int*num_times
        print(case,' num_times: ',num_times)
        time_sec=np.linspace(time_beg,time_end,num_times)
        run_dict['df']=df
        run_dict['L0']=L0
        run_dict['N']=N
        run_dict['df']['time_sec']=time_sec
        run_dict['df']['time_nd']=time_sec*N
        case_list.append((case,L0))
    case_list.sort(key=lambda case: case[1])
    plt.close('all')
    plotlist=['h','h1','h0','delhtop','delhbot','delhtot','delgm',
              'delhtot_rt','zf','zf0','zf1','delzfbot']
    xy_dict={'h':('time_nd','h_nd'),'h1':('time_nd','h1_nd'),'h0':('time_nd','h0_nd'),'delhtop':('time_nd','delhtop'),
             'delhbot':('time_nd','delhbot'),'delhtot':('time_nd','delhtot'),'delgm':('time_nd','delgm'),
             'delhtot_rt':('time_sec','delhtot_rt'),'zf':('time_nd','zf_nd'),
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
        df=run_key[casename]['df']
        zenc=find_zenc(df['time_sec'],N,L0)
        for key in ['h','h0','h1','zf','zf0','zf1']:
            nd_key='{}_nd'.format(key)
            df[nd_key]=df[key]/zenc
        df['delhtop']=(df['h1'] - df['h'])/df['h']
        df['delhbot']=(df['h'] - df['h0'])/df['h']
        df['delhtot']=(df['h1'] - df['h0'])/df['h']
        df['delgm']=0.55*(zenc/L0)**(-2/3.)  #eq. 26
        df['delhtot_rt']=(df['h1'] - df['h0'])/df['h']
        df['delzfbot']=(df['zf'] - df['zf0'])/df['zf']
        df['zenc']=zenc
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
    df_cases['N']=N
    df_cases['fluxes']=fluxes
    df_cases['gammas']=gammas
    for key,axis in plot_dict.items():
        filename='{}.png'.format(key)
        axis.figure.savefig(filename)
    plt.show()
    with pd.HDFStore('paper_table.h5','w') as store:
        store.put('cases',df_cases,format='table')
        store.get_storer('cases').attrs.history = 'written 2015/8/5'
        



#Files = /newtera/tera/phil/nchaparr/python/Plotting/rundate/data/AvProfLims
#indices of h0, h, h1 = 0, 1, 2
#rundate = [Dec142013, Nov302013, Dec202013, Dec252013, Jan152014_1, Mar12014, Mar52014]
#sfcflx/gamma = [100/10, 100/5, 60/5, 60/2.5, 150/5, 60/10, 150/10]
#time increments for all except Nov302013 = 600, for Nov302013 = 900
#end time = 28800

