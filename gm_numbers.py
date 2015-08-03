run_key={'Nov302013':(100,5),         
         'Dec142013':(100,10),        
         'Dec202013':(60,5),          
         'Dec252013':(60,2.5),        
         'Jan152014_1':(150,5),         
         'Mar12014':(60,10),         
         'Mar52014':(150,10)}        

def gm_vars(h,surface_flux,gamma):
    rho=1.
    cp=1004.
    g=9.8
    flux=surface_flux/(rho*cp)  #from W/m^2 to m K/s
    theta_0=280.  #K
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
    return delta,L0

if __name__ == "__main__":
    h=800.
    for case,run_tup in run_key.items():
        surface_flux,gamma=run_tup
        print(case,h,surface_flux,gamma,gm_vars(h,surface_flux,gamma))

#files = /newtera/tera/phil/nchaparr/python/Plotting/rundate/data/AvProfLims
#indices of h0, h, h1 = 0, 1, 2
#rundate = [Dec142013, Nov302013, Dec202013, Dec252013, Jan152014_1, Mar12014, Mar52014]
#sfcflx/gamma = [100/10, 100/5, 60/5, 60/2.5, 150/5, 60/10, 150/10]
#time increments for all except Nov302013 = 600, for Nov302013 = 900
#end time = 28800

