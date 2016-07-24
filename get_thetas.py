"""
  input: json file produced by calc_index.py
  output:  h5 containing all ensemble arrays of wvert and theta perturbations

  python get_thetas.py -j case_details.json -c Dec142013

  which produces thetas_case_name.h5
"""
import numpy as np
import nchap_fun as nc
from Make_Timelist import Make_Timelists
from nchap_class import Get_Var_Arrays1
import h5py
from collections import OrderedDict
import json
from make_tuple import make_tuple
import pdb
from calc_times import list_times
import sys

def get_ensemble(date,dump_time_label):
    """Pull all 10 enssemble theta profiles for all times for a particular date
    """

    filepath = {'root_dir': '/tera/phil/nchaparr/tera2_cp/nchaparr',
                'sam_dir': 'runs/sam_case'}
    suffix = "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_"
    filepath['date'] = date
    prefix = "{root_dir:}/{date:}/{sam_dir:}".format_map(filepath)
    Vars = Get_Var_Arrays1(prefix, suffix, dump_time_label)
    filelist = Vars.get_filelist()
    thetas_list, press_list = Vars.get_thetas(calc_mean=True)
    height = Vars.get_height()
    press = Vars.get_press()
    #filenames = Vars.get_filelist()
    theta_accum = thetas_list[0]
    for count,theta in enumerate(thetas_list[1:]):
        print('theta ensemble: {}'.format(count))
        theta_accum += theta
    theta_accum = theta_accum/len(thetas_list)
    flux_list = Vars.get_wvelthetaperts(calc_mean=True)
    flux_accum = flux_list[0]
    for count,flux in enumerate(flux_list[1:]):
        print('flux ensemble: {}'.format(count))
        flux_accum += flux
    flux_accum += flux
    return filelist, press, height, theta_accum, flux_accum


def write_h5(case_dict,theta_array,flux_array,press,height,h5_outfile='test.h5'):
    with h5py.File(h5_outfile, 'w') as f:
        case = f['/']
        dset = case.create_dataset('thetas',theta_array.shape, dtype=theta_array.dtype)
        dset[...] = theta_array[...]
        dset = case.create_dataset('fluxes',flux_array.shape, dtype=flux_array.dtype)
        dset[...] = flux_array[...]
        dset = case.create_dataset('height',height.shape, dtype=height.dtype)
        dset[...]=height[...]
        dset = case.create_dataset('press',press.shape, dtype=press.dtype)
        dset[...]=press[...]
        case.attrs['case_dict'] = json.dumps(case_dict,indent=4)
    return None


if __name__ == "__main__":

    import argparse, textwrap
    linebreaks = argparse.RawTextHelpFormatter
    descrip = textwrap.dedent(globals()['__doc__'])
    parser = argparse.ArgumentParser(formatter_class=linebreaks,
                                     description=descrip)
    parser.add_argument('-j', '--jfile', help='json file with run info', required=True)
    parser.add_argument('-c',
                        '--case',
                        help='case name',
                        required=True)
    parser.add_argument('-l',
                    '--case_list',
                    help='boolean dump case names',
                        action = 'store_true',
                        required=False)

    args = parser.parse_args()

    with open(args.jfile,'r') as f:
        case_dict = json.load(f)
        if args.case_list:
            print(case_dict.keys())
            sys.exit(1)
    case_list = ['Mar52014', 'Dec252013', 'Dec202013', 'Mar12014', 'Jan152014_1', 'Nov302013']
    for case_name in case_list:    
        time_tup, dump_time_list, Times = list_times(case_name)
        theta_list=[]
        flux_list = []
        for time_stamp in dump_time_list:
            print('timestamp {}'.format(time_stamp))
            filelist, press, height, thetas, fluxes = get_ensemble(case_name,time_stamp)
            theta_list.append(thetas)
            flux_list.append(fluxes)

        theta_array = np.array(theta_list)
        flux_array = np.array(flux_list)
        the_case = case_dict[case_name]
        the_case['time_strings'] = dump_time_list
        the_case['filelist'] = filelist
        the_case['float_hours'] = list(Times)
        outfile ='thetaprofs_{}.h5'.format( case_name)
        write_h5(the_case,theta_array,flux_array,press,height,h5_outfile=outfile)


