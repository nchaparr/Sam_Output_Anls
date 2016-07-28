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
    nensembles = len(thetas_list)
    for count,theta in enumerate(thetas_list[1:]):
        print('theta ensemble: {}'.format(count))
        theta_accum += theta
    theta_avg = theta_accum/nensembles
    meanflux_list, quad_list = Vars.get_wvelthetaperts(calc_mean=True,quadrants=True)
    ensemble_list = list(zip(meanflux_list, quad_list))
    mean_accum = meanflux_list[0]
    quad_accum = quad_list[0]
    for mean_flux, quad_flux in ensemble_list:
        print('flux ensemble: {}'.format(count))
        mean_accum += mean_flux
        quad_accum += quad_flux
    mean_flux_avg = mean_accum/nensembles
    quad_flux_avg = quad_accum/nensembles

    return filelist, press, height, theta_avg, mean_flux_avg, quad_flux_avg


def write_h5(case_dict,theta_array,mean_flux_array,quad_flux_array,press,height,h5_outfile='test.h5'):
    with h5py.File(h5_outfile, 'w') as f:
        case = f['/']
        dset = case.create_dataset('thetas',theta_array.shape, dtype=theta_array.dtype)
        dset[...] = theta_array[...]
        dset = case.create_dataset('mean_fluxes',mean_flux_array.shape, dtype=mean_flux_array.dtype)
        dset[...] = mean_flux_array[...]
        dset = case.create_dataset('quad_fluxes',quad_flux_array.shape, dtype=quad_flux_array.dtype)
        dset[...] = quad_flux_array[...]
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
    case_list = ['Mar12014','Nov302013','Mar52014', 'Dec142013','Dec252013', 'Dec202013',  'Jan152014_1']
    for case_name in case_list[:1]:    
        time_tup, dump_time_list, Times = list_times(case_name)
        theta_list=[]
        mean_flux_list = []
        quad_flux_list = []
        for time_stamp in dump_time_list:
            print('timestamp {}'.format(time_stamp))
            filelist, press, height, thetas, mean_flux, quad_fluxes = get_ensemble(case_name,time_stamp)
            theta_list.append(thetas)
            mean_flux_list.append(mean_flux)
            quad_flux_list.append(quad_fluxes)

        theta_array = np.array(theta_list)
        mean_flux_array = np.array(mean_flux_list)
        quad_flux_array = np.array(quad_flux_list)
        
        the_case = case_dict[case_name]
        the_case['time_strings'] = dump_time_list
        the_case['filelist'] = filelist
        the_case['float_hours'] = list(Times)
        outfile ='thetaprofs_{}.h5'.format( case_name)
        write_h5(the_case,theta_array,mean_flux_array, quad_flux_array,press,height,h5_outfile=outfile)


