import sys
from scipy.integrate import ode
import numpy as np
from datetime import datetime
import scipy.special as scsp 
import time

"""
   Zero order CBL Mixed layer model
   
"""

#Constants, given parameters
cp = 1004 #specific heat capacity J Kg-1 K-1
rho  = 1.225 #density of air Kg m-3
w_h = 0 # subsidence in ma s-1
k = .3 # need to check what this constant is exactly
Flx = 60 #constant surface flux w m-2
gamma = 1.0*5/1000 #lapse rate above the inversion K m-1
theta0 = 301.448318781 
height0 = 580.0
jump0 = 0.696655273438

def MLZero():    
    yinit = [theta0, height0, jump0]
    tinit = 3540
    tfin = 28800
    dt = 60
   
    r = ode(F).set_integrator('vode', method = 'BDF', order=15, nsteps = 300000)
    r.set_f_params()      
    r.set_initial_value(yinit, tinit)
   
    y = np.array(yinit)
    t = np.array(tinit)

    Start = time.clock()
    
    while r.successful() and r.t < tfin:
        r.integrate(r.t+dt)
        print 'now at', r.t
        
        Timediff = time.clock()-Start
        Start = time.clock()        
        dump(r.y, '/tera/phil/nchaparr/python/Plotting/Sep302013/data/MLZero_Vars.txt')
        dump(F(r.t, r.y), '/tera/phil/nchaparr/python/Plotting/Sep302013/data/MLZero_Diffs.txt')
        dump([r.t, Timediff], '/tera/phil/nchaparr/python/Plotting/Sep302013/data/MLZero_Checks.txt')
        
def F(t, y):
    yp = []
    #time_hrs = 1.0*t/3600
    #Flx = scsp.erf(time_hrs/(2.5*np.sqrt(2)))
    Vars = Calc_Vars(k, Flx, cp, rho, y[1],  y[2], w_h, gamma)
    yp.append(Vars[0]) 
    yp.append(Vars[1])
    yp.append(Vars[2])
    return yp

def Calc_Vars(k, Flx, cp, rho, h, jump, w_h, gamma):
    dtheta = 1.0*(1+k)*Flx/(h*cp*rho)
    dheight = 1.0*(k*Flx)/(rho*cp*jump) + w_h
    djump = (dheight - w_h)*gamma - dtheta
    Vars=[dtheta, dheight, djump]
    return Vars

def dump(array, filename):
    a = "  ".join(map(str, array))
    output = open(filename, 'ab')
    output.write("%s\n"%str(a))        
    output.close()

if __name__ == "__main__":
    MLZero()







