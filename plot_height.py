import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from Make_Timelist import *
import sys
sys.path.insert(0, '/tera/phil/nchaparr/python')
import nchap_fun as nc
from matplotlib.lines import Line2D

from matplotlib import rcParams
rcParams.update({'font.size': 10})
"""plots scaled dhdt vs inverse richardson number based on average profile
   now has option to do least squares, and plot 
"""

dump_time_list, Times = Make_Timelists(1, 600, 28800)
Times = np.array(Times)
dump_time_list0, Times0 = Make_Timelists(1, 900, 28800)

dhdtplot = np.genfromtxt(
    "/tera/users/nchaparr/Dec142013/data/dhdtinvriplt.txt")
dhdtplot0 = np.genfromtxt(
    "/tera/users/nchaparr/Nov302013/data/dhdtinvriplt.txt")
dhdtplot1 = np.genfromtxt(
    "/tera/users/nchaparr/Dec202013/data/dhdtinvriplt.txt")
dhdtplot2 = np.genfromtxt(
    "/tera/users/nchaparr/Dec252013/data/dhdtinvriplt.txt")
dhdtplot3 = np.genfromtxt(
    "/tera/users/nchaparr/Jan152014_1/data/dhdtinvriplt.txt")
dhdtplot4 = np.genfromtxt(
    "/tera/users/nchaparr/Mar12014/data/dhdtinvriplt.txt")
dhdtplot5 = np.genfromtxt(
    "/tera/users/nchaparr/Mar52014/data/dhdtinvriplt.txt")

prefix = './fig6c/'
AvProfVars = np.genfromtxt(prefix + 'Dec142013/data/AvProfLims')
AvProfVars0 = np.genfromtxt(prefix + 'Nov302013/data/AvProfLims')
AvProfVars1 = np.genfromtxt(prefix + 'Dec202013/data/AvProfLims')
AvProfVars2 = np.genfromtxt(prefix + 'Dec252013/data/AvProfLims')
AvProfVars3 = np.genfromtxt(prefix + 'Jan152014_1/data/AvProfLims')
AvProfVars4 = np.genfromtxt(prefix + 'Mar12014/data/AvProfLims')
AvProfVars5 = np.genfromtxt(prefix + 'Mar52014/data/AvProfLims')

rinovals = np.genfromtxt(prefix + 'Dec142013/data/invrinos')
rinovals0 = np.genfromtxt(prefix + 'Nov302013/data/invrinos')
rinovals1 = np.genfromtxt(prefix + 'Dec202013/data/invrinos')
rinovals2 = np.genfromtxt(prefix + 'Dec252013/data/invrinos')
rinovals3 = np.genfromtxt(prefix + 'Jan152014_1/data/invrinos')
rinovals4 = np.genfromtxt(prefix + 'Mar12014/data/invrinos')
rinovals5 = np.genfromtxt(prefix + 'Mar52014/data/invrinos')

Deltah0 = np.subtract(AvProfVars[:, 2], AvProfVars[:, 0])
Deltah0 = np.divide(Deltah0, AvProfVars[:, 1])
Deltah00 = np.subtract(AvProfVars0[:, 2], AvProfVars0[:, 0])
Deltah00 = np.divide(Deltah00, AvProfVars0[:, 1])
Deltah01 = np.subtract(AvProfVars1[:, 2], AvProfVars1[:, 0])
Deltah01 = np.divide(Deltah01, AvProfVars1[:, 1])
Deltah02 = np.subtract(AvProfVars2[:, 2], AvProfVars2[:, 0])
Deltah02 = np.divide(Deltah02, AvProfVars2[:, 1])
Deltah03 = np.subtract(AvProfVars3[:, 2], AvProfVars3[:, 0])
Deltah03 = np.divide(Deltah03, AvProfVars3[:, 1])
Deltah03[15:21] = np.nan
Deltah04 = np.subtract(AvProfVars4[:, 2], AvProfVars4[:, 0])
Deltah04 = np.divide(Deltah04, AvProfVars4[:, 1])
Deltah04[10:17] = np.nan
Deltah05 = np.subtract(AvProfVars5[:, 2], AvProfVars5[:, 0])
Deltah05 = np.divide(Deltah05, AvProfVars5[:, 1])

#plot the heights vs time
#print Line2D.markers
Fig2 = plt.figure(2)
Fig2.clf()
Ax3 = Fig2.add_subplot(111)
FitFunc = np.polyfit(Times[11:], AvProfVars5[11:, 1], 2, full=False)
Fit = FitFunc[0] * Times[11:]**2 + FitFunc[1] * Times[11:] + FitFunc[2]
dhdt = 1.0 * (2 * FitFunc[0] * Times[11:] + FitFunc[1]) / 3600
deltah = np.subtract(AvProfVars5[:, 2], AvProfVars5[:, 0])
deltah = np.divide(deltah, AvProfVars5[:, 1])
tau = 1.0 * rinovals5[:, 4] / 3600
scaled_time = np.divide(Times, tau)
scaled_dhdt = np.divide(dhdt, rinovals5[11:, 2])
dhdtinvriplt = np.vstack((rinovals5[11:, 1], scaled_dhdt))
dhdtinvriplt = np.transpose(np.vstack((dhdtinvriplt, deltah[11:])))


Ax3.loglog(rinovals[11:, 1], Deltah0[11:], 'kv', label='100/10', markersize=10)
Ax3.loglog(rinovals0[7:, 1], Deltah00[7:], 'ko', label='100/5', markersize=10)
Ax3.loglog(rinovals1[11:, 1], Deltah01[11:], 'yo', label='60/5', markersize=10)
Ax3.loglog(rinovals2[11:, 1],
           Deltah02[11:],
           'y*',
           label='60/2.5',
           markersize=10)
Ax3.loglog(rinovals3[11:29, 1],
           Deltah03[11:29],
           'ro',
           label='150/5',
           markersize=10)
Ax3.loglog(rinovals4[11:, 1],
           Deltah04[11:],
           'yv',
           label='60/10',
           markersize=10)
Ax3.loglog(rinovals5[11:, 1],
           Deltah05[11:],
           'rv',
           label='150/10',
           markersize=10)

Ax3.set_xlabel(r"$Ri_{\Delta}^{-1}$", fontsize=30)
Ax3.tick_params(axis="both", labelsize=18)
Ax3.text(.086, .87, r'(c)', fontsize=30)
#Ax3.set_ylabel(r"$h \ (m)$", fontsize=20)
plt.ylim(.2, 1)
plt.xlim(.02, .1)
Ax3.set_xticks([.02, .04, .06, .08, .1])
Ax3.set_xticklabels([0.02, 0.04, 0.06, 0.08, 0.1])
Ax3.set_yticks([.2, .4, .6, .8, 1])
Ax3.set_yticklabels([0.2, 0.4, 0.6, 0.8, 1])
plt.tight_layout()
plt.show()
