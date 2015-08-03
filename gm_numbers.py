def gm_vars(h,surface_flux,gamma):
    rho=1.
    cp=1004.
    g=9.8
    surface_flux=60 #W/m^2
    flux=surface_flux/(rho*cp)  #m K/s
    theta_0=280.  #K
    B0=flux*g/theta_0
    gamma=2.5/1000.  #K/m
    g=9.8  #m/s^2
    N2=g/theta_0*gamma  #s**(-2)
    N=N2**0.5
    L0=(B0/N**3.)**0.5  #gm eqn 3
    Re0=(L0*B0)**(1./3.)
    h=800.
    wstar=(g*h/theta_0*flux)**(1./3)
    thetastar=flux/wstar
    wstar_gm=(B0*h)**(1./3.)
    c_gamma=0.55
    delta=c_gamma*wstar/N

if __name__ == "__main__":



