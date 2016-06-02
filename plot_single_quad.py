"""
plot individual quadrant profiles for times from 300-500
and compare to the total flux at those times

example:  python plot_single_quad.py -q flux_out_300.h5

"""

import h5py
from matplotlib import pyplot as plt
from collections import defaultdict
import json
import numpy as np
from gm_numbers import find_zenc
import glob


if __name__ == "__main__":

    import argparse, textwrap
    linebreaks = argparse.RawTextHelpFormatter
    descrip = textwrap.dedent(globals()['__doc__'])
    parser = argparse.ArgumentParser(formatter_class=linebreaks,
                                     description=descrip)
    flux_files= glob.glob('flux*h5')
    

    attr_dict = dict(hvals='height_columns',scales='scale_columns')
    big_dict={}
    for flux_file in flux_files:
        with h5py.File(flux_file,'r') as infile:
            print('opening {}'.format(flux_file))
            case_numbers = json.loads(infile['/'].attrs['case_dict'])
            _,_,t_nd = flux_file.split('_')
            t_nd,_ = t_nd.split('.')
            t_nd = int(t_nd)
            print(case_numbers)
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
        big_dict[t_nd] = data_dict
        
    plot = True
    case_keys =  ['Mar12014','Dec142013','Mar52014']
    quad_keys = ['wn_tp', 'wn_tn','wp_tp','wp_tn']
    if plot:
        plt.close('all')
        for case in case_keys:
            for quad in quad_keys:
                fig,ax = plt.subplots(1,1)
                times = list(big_dict.keys())
                times.sort()
                for t_nd in times:
                    data_dict = big_dict[t_nd]
                    print('debug 0',data_dict.keys())
                    height_columns = data_dict[case]['hvals']['height_columns']
                    height = data_dict[case]['height']['data'][...]
                    zg0_index = data_dict[case]['hvals']['height_columns']['zg0']
                    zf0_index = data_dict[case]['hvals']['height_columns']['zf0']
                    time_index = data_dict[case]['time_index']
                    hvals = data_dict[case]['hvals']['data']
                    time_sec = data_dict[case]['time_seconds']
                    N = case_numbers[case]['N']
                    L0 = case_numbers[case]['L0']
                    print('debug: ',time_sec,N,L0)
                    zenc = find_zenc(time_sec,N,L0)
                    height_nd = height/zenc
                    print('found zenc: ',zenc)
                    total_flux=np.zeros_like(height)
                    for key in quad_keys:
                        flux = data_dict[case][key]['data']
                        total_flux = total_flux + flux
                    quad_flux = data_dict[case][quad]['data']
                    ax.plot(quad_flux,height_nd,label=t_nd)
                    ax.plot(total_flux,height_nd)
                    ax.plot(total_flux,height_nd,label='total_flux')
                    ax.axhline(hvals[time_index,zg0_index]/zenc)
                    ax.axhline(hvals[time_index,zf0_index]/zenc)
                title = "quadrant: {}, case: {}".format(quad,case)
                ax.set(title=title,ylim=(0.,1.5))
                figname = 'timeseries_{}_{}.png'.format(case,quad)
                ax.legend()
                fig.savefig(figname)

        plt.show()
