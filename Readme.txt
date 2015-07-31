Run Key

Surface Heat Flux/Gamma            Run Date
100/5                              Nov302013
100/10                             Dec142013
60/5                               Dec202013
60/2.5                             Dec252013
150/5                              Jan152014_1
60/10                              Mar12014
150/10                             Mar52014                         

Most recent sam code (test on high resolution the whole way up) is in 

users/nchaparr/sam_grex

the grd file needs to be edited for each run, higher vertical region of high resolution for lower gamma/higher sfc heatflux

Code for editing .nc input files is in

newtera/tera/phil/nchaparr/python/Pert_Files
1. make_snd.py
2. snd2nc1.py
3. snd2nc2.py

ensemble code and output are in

/newtera/tera/phil/nchaparr/sam_grex_ensemble (backed up in /newtera/tera/phil/nchaparr/Aug122014)
