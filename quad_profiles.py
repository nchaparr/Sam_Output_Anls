"""
  Given a set of profiles produces by Flux_Quads.py  oputput conditionally sampled vertical
  quadrant profiles

  Example:  python quad_profiles.py -p profiles_may26_300.h5 -o flux_out_300.h5

"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
from  collections import defaultdict
from enum import IntEnum

class Heights(IntEnum):
    zg0 = 0
    zf0 = 3


runs = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

if __name__ == "__main__":

    import argparse, textwrap
    linebreaks = argparse.RawTextHelpFormatter
    descrip = textwrap.dedent(globals()['__doc__'])
    parser = argparse.ArgumentParser(formatter_class=linebreaks,
                                     description=descrip)
    parser.add_argument('-p', '--prof', help='h5 file produced by Flux_Quads.py', required=True)
    parser.add_argument('-o',
                        '--out',
                        help='h5 file containing the quadrant output',
                        required=True)
    args = parser.parse_args()

    
    h5_profiles = args.prof
    flux_profiles = args.out
    case_dict={}
    
    keys=['wntp','wntn','wptp','wptn','wt','dwtdz', 'wn_tp','wn_tn','wp_tn','wp_tp','tp_wn', 'tp_wp','tn_wn','tn_wp','dwn_tpdz','dwn_tndz','dwp_tndz','dwp_tpdz','dtp_wndz','dtp_wpdz','dtn_wndz','dtn_wpdz']
    keys=keys[:4]
    zeros = np.zeros([len(keys)]) #chang to 5
    
    with h5py.File(h5_profiles,'r') as infile, h5py.File(flux_profiles,'w') as outfile:
        firstpass = True
        for case in list(infile.keys()):
            case_group = outfile.create_group(case)
            case_dict[case] = {}
            case_dict[case]['flux'] = []
            
            height = infile[case]['height'][...]
            
            dheight = np.diff(height)
            element0=[0]
            dheight = np.hstack((element0, dheight))
            
            
            hvals = infile[case]['hvals'][...]
            print(infile[case]['time_list'].attrs['time_index'])
            time_index = int(infile[case]['time_list'].attrs['time_index'])
            zg0 = hvals[time_index,int(Heights.zg0)]
            zf0 = hvals[time_index,int(Heights.zf0)]
            case_dict[case]['zg0'] = zg0
            case_dict[case]['zf0'] = zf0
            case_dict[case]['height'] = height
            for key in keys:
                case_dict[case][key]=[]
            flux_prof = []
            print(case)
            if firstpass:
                #
                # reuse these arrays so memory doesn't grow 
		
                nz,ny,nx = infile[case]['0']['thetapert'].shape
                thetapertslice = np.empty_like(infile[case]['0']['thetapert'][0,...])
                wvelpertslice = np.empty_like(infile[case]['0']['thetapert'][0,...])
                wvelthetapertslice = np.empty_like(infile[case]['0']['thetapert'][0,...])
		              	        
                firstpass = False
            #for run in runs:
            for lev in range(nz):  
                #print(keys)
                store_sum = dict(zip(keys,zeros))
                #print(store_sum)

                for run in runs:
                    thetapert=np.empty_like(infile[case][run]['thetapert'][...])
                    thetapert[...] = infile[case][run]['thetapert'][...]                    
                    dthetapert=np.diff(thetapert,axis=0)                    
                    element0 = np.zeros_like(thetapert[0:1,:,:])                    
                    dthetapert=np.concatenate((element0, dthetapert),axis=0)
                    
                
                    wvelpert=np.empty_like(infile[case][run]['wvelpert'][...])
                    wvelpert[...] = infile[case][run]['wvelpert'][...]
                    dwvelpert=np.diff(wvelpert,axis=0)
                    dwvelpert=np.concatenate((element0, dwvelpert), axis=0) 
                    
                    wvelthetapert=np.empty_like(infile[case][run]['wvelthetapert'][...])
                    wvelthetapert[...] = infile[case][run]['wvelthetapert'][...]
                    dwvelthetapert=np.diff(wvelthetapert,axis=0)
                    dwvelthetapert=np.concatenate((element0, dwvelthetapert), axis=0) 
            
                
                    thetapertslice[...] = infile[case][run]['thetapert'][lev,...]
                    wvelpertslice[...] = infile[case][run]['wvelpert'][lev,...]
                    wvelthetapertslice[...]=infile[case][run]['wvelthetapert'][lev,...] 
                    dthetapertslice=dthetapert[lev,:,:]
                    dwvelpertslice=dwvelpert[lev,:,:]
                    dwvelthetapertslice=dwvelthetapert[lev,:,:]
                    
                    hit = np.logical_and(wvelpertslice < 0.,thetapertslice > 0.)
                    flux = (thetapertslice[hit]).mean()*(wvelpertslice[hit]).mean()
                    store_sum['wntp'] += flux
                    flux = np.sqrt((wvelpertslice[hit]**2).mean())
                    store_sum['wn_tp'] += flux
                    flux=(dwvelpertslice[hit]).mean()
                    store_sum['dwn_tpdz']+=flux
                    flux = np.sqrt((thetapertslice[hit]**2).mean())
                    store_sum['tp_wn'] += flux
                    flux=(dthetapertslice[hit]).mean()
                    store_sum['dtp_wndz']+=flux


                    hit = np.logical_and(wvelpertslice < 0.,thetapertslice < 0.)
                    flux = (thetapertslice[hit]).mean()*(wvelpertslice[hit]).mean()
                    store_sum['wntn'] += flux
                    flux = np.sqrt((wvelpertslice[hit]**2).mean())
                    store_sum['wn_tn'] += flux
                    flux=(dwvelpertslice[hit]).mean()
                    store_sum['dwn_tndz'] += flux  
                    flux = np.sqrt((thetapertslice[hit]**2).mean())
                    store_sum['tn_wn'] += flux
                    flux=(dthetapertslice[hit]).mean()
                    store_sum['dtn_wndz'] += flux  


                    hit = np.logical_and(wvelpertslice > 0.,thetapertslice < 0.)
                    flux = (thetapertslice[hit]).mean()*(wvelpertslice[hit]).mean()
                    store_sum['wptn'] += flux                    
                    flux = np.sqrt((wvelpertslice[hit]**2).mean())
                    store_sum['wp_tn'] += flux
                    flux=(dwvelpertslice[hit]).mean()
                    store_sum['dwp_tndz'] += flux  
                    flux = np.sqrt((thetapertslice[hit]**2).mean())
                    store_sum['tn_wp'] += flux
                    flux=(dthetapertslice[hit]).mean()
                    store_sum['dtn_wpdz'] += flux  


                    hit = np.logical_and(wvelpertslice > 0.,thetapertslice > 0.)
                    flux = (thetapertslice[hit]).mean()*(wvelpertslice[hit]).mean()
                    store_sum['wptp'] += flux
                    flux = np.sqrt((wvelpertslice[hit]**2).mean())
                    store_sum['wp_tp'] += flux
                    flux=(dwvelpertslice[hit]).mean()
                    store_sum['dwp_tpdz']+=flux 
                    flux = np.sqrt((thetapertslice[hit]**2).mean())
                    store_sum['tp_wp'] += flux
                    flux=(dthetapertslice[hit]).mean()
                    store_sum['dtp_wpdz']+=flux 
                    
                    flux=wvelthetapertslice.mean()
                    store_sum['wt'] += flux #
                    flux=dwvelthetapertslice.mean()
                    store_sum['dwtdz'] += flux #


                for key in keys:
                    store_sum[key] = store_sum[key]/len(runs)
                    case_dict[case][key].append(store_sum[key])

            for key in keys:
                flux_prof = np.array(case_dict[case][key]) #expand to include wperts and theta perts
                dset = case_group.create_dataset(key,flux_prof.shape,dtype=flux_prof.dtype)
                dset[:] = flux_prof[:]
            for dname in ['height','hvals','scales']:
                the_var = infile[case][dname][...]
                dset = outfile[case].create_dataset(dname,the_var.shape,the_var.dtype)
                dset[...] = the_var[...]
                for key, value in infile[case][dname].attrs.items():
                    dset.attrs[key] = value
            outfile[case].attrs['time_index'] = time_index
            outfile[case].attrs['time_list'] = infile[case]['time_list'][...]
        outfile['/'].attrs['case_dict'] = infile['/'].attrs['case_dict']
