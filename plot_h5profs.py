import h5py
from matplotlib import pyplot as plt
from collections import defaultdict
import json
import numpy as np
import pandas as pd
from gm_numbers import find_zenc


if __name__ == "__main__":
    paper_table_file = './data/paper_table.h5'
    with pd.HDFStore(paper_table_file,'r') as store:
        df_table = store['cases']
    
    flux_profiles = 'flux_out_300.h5'

    attr_dict = dict(hvals='height_columns',scales='scale_columns')
    with h5py.File(flux_profiles,'r') as infile:
        case_numbers = json.loads(infile['/'].attrs['case_dict'])
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

    keys=['wn_tp', 'wn_tn','wp_tp','wp_tn']
    plot = True
    if plot:
        plt.close('all')
        for case in cases:
            height_columns = data_dict[case]['hvals']['height_columns']
            height = data_dict[case]['height']['data'][...]
            zg0_index = data_dict[case]['hvals']['height_columns']['zg0']
            zf0_index = data_dict[case]['hvals']['height_columns']['zf0']
            time_index = data_dict[case]['time_index']
            hvals = data_dict[case]['hvals']['data']
            time_sec = data_dict[case]['time_seconds']
            N = case_numbers[case]['N']
            L0 = case_numbers[case]['L0']
            zenc = find_zenc(time_sec,N,L0)
            height_nd = height/zenc
            print('found zenc: ',zenc)
            fig,ax = plt.subplots(1,1)
            total_flux=np.zeros_like(height)
            for key in keys:
                flux = data_dict[case][key]['data']
                ax.plot(flux,height_nd,label=key)
                total_flux = total_flux + flux
            ax.plot(total_flux,height_nd,label='total_flux')
            ax.axhline(hvals[time_index,zg0_index]/zenc)
            ax.axhline(hvals[time_index,zf0_index]/zenc)
            title = "time step {time_periods:5.2f} (periods), L0={L0:6.3f} m, Nperiod = {Nperiod:4.2f} min".format_map(case_numbers[case])
            ax.set(title=title,ylim=(0.,1.5))
            figname = '{}_300.png'.format(case)
            ax.legend()
            fig.savefig(figname)

        plt.show()
