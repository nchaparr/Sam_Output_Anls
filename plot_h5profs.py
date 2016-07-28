"""
 Given a set of conditionally sampled profiles produced by quad_profiles.py
produce a profile plot for each case as a png file

example:  python plot_h5profs.py -q flux_out_300.h5

"""

import h5py
from matplotlib import pyplot as plt
from collections import defaultdict
import json
import numpy as np
from gm_numbers import find_zenc


if __name__ == "__main__":

    import argparse, textwrap
    linebreaks = argparse.RawTextHelpFormatter
    descrip = textwrap.dedent(globals()['__doc__'])
    parser = argparse.ArgumentParser(formatter_class=linebreaks,
                                     description=descrip)
    parser.add_argument('-q', '--quad', help='h5 file produced by quad_quad_profiles.py', required=True)
    args = parser.parse_args()
    
    flux_profiles = args.quad

    attr_dict = dict(hvals='height_columns',scales='scale_columns')
    with h5py.File(flux_profiles,'r') as infile:
        case_numbers = json.loads(infile['/'].attrs['case_dict'])
        l=lambda:defaultdict(l)
        data_dict = l()
        cases = list(infile.keys())
        for case in cases:
            group = infile[case]
            data_dict[case]['time_index'] = infile[case].attrs['time_index']
            data_dict[case]['time_list'] =  infile[case].attrs['time_list'][:]
            data_dict[case]['time_seconds'] = data_dict[case]['time_list'][data_dict[case]['time_index']]*3600.
            case_numbers[case]['time_periods'] = data_dict[case]['time_seconds']*case_numbers[case]['N']
            for dset in group.keys():
                data_dict[case][dset]['data']=infile[case][dset][...]
                if dset in ['hvals','scales']:
                    attrname = attr_dict[dset]
                    data_dict[case][dset][attrname] = json.loads(infile[case][dset].attrs[attrname])

    
    #all_keys=['wntp','wntn','wptp','wptn','wt','wn_tp','wn_tn','wp_tn','wp_tp','tp_wn', 'tp_wp','tn_wn','tn_wp','dwn_tpdz','dwn_tndz','dwp_tndz','dwp_tpdz','dtp_wndz','dtp_wpdz','dtn_wndz','dtn_wpdz']
    keys=['wt']
    plot = True
    if plot:
        plt.close('all')
        for key in keys:
            fig,ax = plt.subplots(1,1)

            keysforzs=['max','avg','min']
            zerosforzs=[0,0,0]
            zgdist=dict(zip(keysforzs, zerosforzs))
            zf0dist=dict(zip(keysforzs, zerosforzs))
            zg0dist=dict(zip(keysforzs, zerosforzs))

            for case in cases:
                #total_flux=np.zeros_like(height)
                height_columns = data_dict[case]['hvals']['height_columns']
                height = data_dict[case]['height']['data'][...]                
                time_index = data_dict[case]['time_index']
                scales = data_dict[case]['scales']['data']
                wstar=scales[time_index, 2]
                
                hvals = data_dict[case]['hvals']['data']
                zg=hvals[time_index, 1]
                zg0=hvals[time_index,0]
                zf0=hvals[time_index,3]
                
                if zg>=zgdist['max']:
                    zgdist['max']=zg
                
                zgdist['avg']+=zg/len(cases)
                
                if zg<=zgdist['min']:
                    zgdist['min']=zg


                if zg0>=zg0dist['max']:
                    zg0dist['max']=zg0
                
                zg0dist['avg']+=zg0/len(cases)
                
                if zg0<=zg0dist['min']:
                    zg0dist['min']=zg0

                if zf0>=zf0dist['max']:
                    zf0dist['max']=zf0
                
                zf0dist['avg']+=zf0/len(cases)
                
                if zf0<=zf0dist['min']:
                    zf0dist['min']=zf0   
    

                thetastar=scales[time_index, 9]#[rino, invrino, wstar, S, tau, mltheta, deltatheta, pi3, pi4, thetastar, c_delta]
                
                time_sec = data_dict[case]['time_seconds']
                N = case_numbers[case]['N']
                L0 = case_numbers[case]['L0']                
                surface_flux=case_numbers[case]['fluxes']/(1004*1.14)
                zenc = find_zenc(time_sec,N,L0)
                print(time_sec, time_index, zg, zg/zenc)
                height_nd = height/zg
                legend=case_numbers[case]['legends']
                #print('found zenc: ',zenc, scales[time_index, 2])	    
                flux = data_dict[case][key]['data']
                #print(flux.shape)
                ax.plot(1.0*flux/surface_flux,height_nd,legend, markersize=10, label=int(L0))
                #ax.axhline(hvals[time_index,zg0_index]/zg)
                
            title = key
            plt.xlabel("scaled total heat flux", size=20)
            plt.ylabel("scaled height",size=20)
            ax.set(title="",ylim=(0,1.2), xlim=(-.2, 1))
            #figname = '{}_100.png'.format(key)
            #ax.legend(numpoints=1, loc='best')
            

            ax.plot([-.2,.7], [zgdist['avg']/zgdist['avg'], zgdist['avg']/zgdist['avg']], 'k:')
            ax.plot([-.2,.7], [zg0dist['avg']/zgdist['avg'], zg0dist['avg']/zgdist['avg']], 'k:')
            ax.plot([-.2,.7], [zf0dist['avg']/zgdist['avg'], zf0dist['avg']/zgdist['avg']], 'k:')
            
            ax.text(.8, zgdist['avg']/zgdist['avg'], r"$\overline{z_{g}}$", size=30)
            ax.text(.8, zf0dist['avg']/zgdist['avg'], r"$\overline{z_{f0}}$", size=30)
            ax.text(.8, zg0dist['avg']/zgdist['avg'], r"$\overline{z_{g0}}$", size=30)
            
            #ax.axhline(zgdist['avg']/zgdist['avg'])
            #ax.axhline(zg0dist['avg']/zgdist['avg'])
            #ax.axhline(zf0dist['avg']/zgdist['avg'])
            #fig.savefig(figname)
                
               
            
        plt.show()
