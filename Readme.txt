Run Key

Surface Heat Flux/Gamma            Run Date
100/5                              Nov302013
100/10                             Dec142013
60/5                               Dec202013
60/2.5                             Dec252013
150/5                              Jan152014_1
60/10                              Mar12014
150/10                             Mar52014                         

__________________________________________________________________________________________________________________________-

Most recent sam code (test on high resolution the whole way up) is in 

users/nchaparr/sam_grex

the grd file needs to be edited for each run, higher vertical region of 
high resolution for lower gamma/higher sfc heatflux

______________________________________________________________________________________________________________________________

Code for editing .nc input files is in

newtera/tera/phil/nchaparr/python/Pert_Files
1. make_snd.py
2. snd2nc1.py
3. snd2nc2.py

____________________________________________________________________________________________________________________

ensemble code and output are in

/newtera/tera/phil/nchaparr/sam_grex_ensemble 
(backed up in /newtera/tera/phil/nchaparr/tera2_cp/nchaparr/Aug122014)

____________________________________________________________________________________________________________________

data txt files: 

output > Sam_Output_Anls > /newtera/tera/phil/nchaparr/python/Plotting/Run Date/data/

eg

	get_limits.py > AvProfLims = [h0, h, h1, zf0, zf, zf1, deltatheta (based on av theta profile), mltheta (av mixed layer pot temp based on av theta profile)]), invrinos [nchap_fun.py: calc_rino(BLHeight, MLTheta, SfcFlux, Theta_jump, gamma, delta_h)], 
	
	Ens_Profs.py > heights0000028800, theta_bar0000028800, wvelthetapert0000000600..
        
        Flux_Quads.py > flux_quads0000001800...       
        
	
	time increments for all except Nov302013 = 600, for Nov302013 = 900
	end time = 28800

_______________________________________________________________________________________________________________________	
