"""
  Given a set of profiles produces by Flux_Quads.py  oputput conditionally sampled vertical
  quadrant profiles

  Example:  python full_fluxes.py -p profs_with_avgtheta_300.h5 -o full_out_300.h5

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

    outvars=['wntp','wntn','wptp','wptn','wt','dwtdz', 'wn_tp','wn_tn','wp_tn','wp_tp','tp_wn', 'tp_wp','tn_wn','tn_wp','dwn_tpdz','dwn_tndz','dwp_tndz','dwp_tpdz','dtp_wndz','dtp_wpdz','dtn_wndz','dtn_wpdz']

    
    with h5py.File(h5_profiles,'r') as infile, h5py.File(flux_profiles,'w') as outfile:
        h5_levs=defaultdict(dict)
        h5_cases=dict()
        nz=312
        case_list=list(infile.keys())[:1]
        outvars=outvars[:4]
        outvars.extend(['wpert','thetapert'])
        height_levs=list(range(nz))
        for case in case_list:
            case_group = outfile.create_group(case)
            h5_cases[case]=case_group
            for lev in height_levs:
                lev_group=case_group.create_group(str(lev))
                h5_levs[case][lev]=lev_group
        firstpass = True
        for case in case_list:
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
            for key in outvars:
                case_dict[case][key]=[]
            flux_prof = []
            print(case)
            if firstpass:
                #
                # reuse these arrays so memory doesn't grow 
                
                nz,ny,nx = infile[case]['0']['thetapert'].shape
                thetapertslice = np.empty_like(infile[case]['0']['thetapert'][0,...]).ravel()
                wvelpertslice = np.empty_like(infile[case]['0']['thetapert'][0,...]).ravel()
                wvelthetapertslice = np.empty_like(infile[case]['0']['thetapert'][0,...]).ravel()
                firstpass = False
            #for run in runs:
            for lev in height_levs:
                print('height lev is: ',lev)
                #print(outvars)
                store_sum = defaultdict(list)


                for run in runs:
                    print(lev,height[lev],run)
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
            
                
                    thetapertslice[...] = infile[case][run]['thetapert'][lev,...].ravel()
                    wvelpertslice[...] = infile[case][run]['wvelpert'][lev,...].ravel()
                    wvelthetapertslice[...]=infile[case][run]['wvelthetapert'][lev,...].ravel()
                    print('initial array size is: ',thetapertslice.shape)
                    dthetapertslice=dthetapert[lev,:,:]
                    dwvelpertslice=dwvelpert[lev,:,:]
                    dwvelthetapertslice=dwvelthetapert[lev,:,:]
                    
                    hit = np.logical_and(wvelpertslice < 0.,thetapertslice > 0.)
                    flux = (thetapertslice[hit])*(wvelpertslice[hit])
                    store_sum['wntp'].append(flux)
                    # flux = np.sqrt(wvelpertslice[hit]**2)
                    # store_sum['wn_tp'].append(flux)
                    # flux=(dwvelpertslice[hit])
                    # store_sum['dwn_tpdz'].append(flux)
                    # flux = np.sqrt(thetapertslice[hit]**2)
                    # store_sum['tp_wn'].append(flux)
                    # flux=(dthetapertslice[hit])
                    # store_sum['dtp_wndz'].append(flux)


                    hit = np.logical_and(wvelpertslice < 0.,thetapertslice < 0.)
                    flux = (thetapertslice[hit])*(wvelpertslice[hit])
                    store_sum['wntn'].append(flux)
                    # flux = np.sqrt(wvelpertslice[hit]**2)
                    # store_sum['wn_tn'].append(flux)
                    # flux=(dwvelpertslice[hit])
                    # store_sum['dwn_tndz'].append(flux)  
                    # flux = np.sqrt((thetapertslice[hit]**2))
                    # store_sum['tn_wn'].append(flux)
                    # flux=(dthetapertslice[hit])
                    # store_sum['dtn_wndz'].append(flux)  


                    hit = np.logical_and(wvelpertslice > 0.,thetapertslice < 0.)
                    flux = (thetapertslice[hit])*(wvelpertslice[hit])
                    store_sum['wptn'].append(flux)                    
                    # flux = np.sqrt(wvelpertslice[hit]**2)
                    # store_sum['wp_tn'].append(flux)
                    # flux=(dwvelpertslice[hit])
                    # store_sum['dwp_tndz'].append(flux)  
                    # flux = np.sqrt(thetapertslice[hit]**2)
                    # store_sum['tn_wp'].append(flux)
                    # flux=(dthetapertslice[hit])
                    # store_sum['dtn_wpdz'].append(flux)  


                    hit = np.logical_and(wvelpertslice > 0.,thetapertslice > 0.)
                    flux = (thetapertslice[hit])*(wvelpertslice[hit])
                    store_sum['wptp'].append(flux)
                    #flux = np.sqrt(wvelpertslice[hit]**2)
                    # store_sum['wp_tp'].append(flux)
                    # flux=(dwvelpertslice[hit])
                    # store_sum['dwp_tpdz'].append(flux) 
                    # flux = np.sqrt(thetapertslice[hit]**2)
                    # store_sum['tp_wp'].append(flux)
                    # flux=(dthetapertslice[hit])
                    # store_sum['dtp_wpdz'].append(flux) 
                    
                    # flux=wvelthetapertslice
                    # store_sum['wt'].append(flux) #
                    # flux=dwvelthetapertslice
                    # store_sum['dwtdz'].append(flux) #
                    store_sum['wpert'].append(wvelpertslice)
                    store_sum['thetapert'].append(thetapertslice)
                #
                # store_sum contains a list of 10 ensemble members at level lev for case case for
                # each key at level l
                #
                for key in outvars:
                    all_points=store_sum[key][0]
                    print('writing: ',case,lev,key,all_points.shape)
                    for ensemble in np.arange(1,10):
                        print(store_sum[key][ensemble].shape)
                        all_points=np.append(all_points,store_sum[key][ensemble])
                    dset = h5_levs[case][lev].create_dataset(key,all_points.shape,dtype=all_points.dtype)
                    dset[...] = all_points[...]
                    
            for dname in ['height','hvals','scales']:
                the_var = infile[case][dname][...]
                dset = outfile[case].create_dataset(dname,the_var.shape,the_var.dtype)
                dset[...] = the_var[...]
                for key, value in infile[case][dname].attrs.items():
                    dset.attrs[key] = value
            outfile[case].attrs['time_index'] = time_index
            outfile[case].attrs['time_list'] = infile[case]['time_list'][...]
        outfile['/'].attrs['case_dict'] = infile['/'].attrs['case_dict']
