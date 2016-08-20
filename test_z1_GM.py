import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)

C=-12822.0429666
gamma=.01
slope1=44.8385331144
thetas = np.arange(0, 400, 10) 
h1 = slope1*thetas + C
h2 = (thetas-300)*1.0/gamma
thetabar = np.genfromtxt("/tera/users/nchaparr/Mar52014/data/" + "theta_bar0000028800")
heights=np.genfromtxt("/tera/users/nchaparr/Mar52014/data/" + "heights0000028800")
theta_0 = heights*gamma + 300



Ax3.plot(thetas, h1, thetas, h2) #theta_0, heights, 
#Ax3.plot([0, 310.103515625, 312],[slope1*(0-310.103515625)+1075, 1075, slope1*(312-310.103515625)+1075])
Ax3.plot(thetabar, heights)
plt.ylim(0, 1200)
plt.xlim(300, 330)
plt.show()
