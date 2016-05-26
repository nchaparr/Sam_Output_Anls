"""
  input: json file produced by calc_index.py
  output:  h5 containing all ensemble arrays of wvert and theta perturbations

  example: python Flux_Quads.py -j index_list_time_nd_300.json  -r may26_300

  which produces profiles_may26_300.h5
"""
import numpy as np
import matplotlib.pyplot as plt
import nchap_fun as nc
from matplotlib.colors import Normalize
from Make_Timelist import Make_Timelists
from nchap_class import Get_Var_Arrays1
import h5py
from collections import OrderedDict
import json
from make_tuple import make_tuple


def get_ensemble(date, dump_time_label, height_level=0):
    """Pulls output from an ensemble cases, gets ensemble averages and perturbations

    Arguments:
    date -- run date eg "Nov302013"
    dump_time_label -- time of output eg '0000000720'
    height_level -- height at which to take slice of perturbations, eg z_g0

    Returns:
    var_bar -- 64 array of horizontally averaged, ensemble averages or perturbations (covariances)
    
    """

    #create lists for variable arrays from each case
    wvelperts_list = []
    thetaperts_list = []

    #get velocity perts and thetas using method from nchap_class
    #sample file name for a case1
    #/tera/phil/nchaparr/tera2_cp/nchaparr/Aug122014/runs/sam_case1/OUT_3D/NCHAPP1_testing_doscamiopdata_24_0000000060.nc
    #
    #hvals order:
    #defined in nchap_fun.Get_CBLHeight

    filepath = {'root_dir': '/tera/phil/nchaparr/tera2_cp/nchaparr',
                'sam_dir': 'runs/sam_case'}
    suffix = "/OUT_3D/keep/NCHAPP1_testing_doscamiopdata_24_"
    filepath['date'] = date
    prefix = "{root_dir:}/{date:}/{sam_dir:}".format_map(filepath)
    Vars = Get_Var_Arrays1(prefix, suffix, dump_time_label)
    file_list = Vars.get_filelist()
    print(file_list)
    thetas_list, press_list = Vars.get_thetas()

    wvels_list = Vars.get_wvelperts()

    height = Vars.get_height()

    filenames = Vars.get_filelist()
    #get arrays of ensemble averaged variables, from nchap_fun

    ens_avthetas = nc.Ensemble1_Average(thetas_list)
    ens_press = nc.Ensemble1_Average(press_list)

    #now get the perturbations
    wvelthetaperts_list = []
    full_profs_list = []
    for i in range(len(
            wvels_list)):  #TODO: this should be more modular, see nchap_class
        thetapert_rough = np.subtract(thetas_list[i], ens_avthetas)
        thetapert = np.zeros_like(thetapert_rough)
        [znum, ynum, xnum] = wvels_list[i].shape
        for j in range(
                znum):  #something like this is done in statistics.f90, staggered grid!
            if j == 0:
                thetapert[j, :, :] = thetapert_rough[j, :, :]
            else:
                thetapert[j, :, :] = 0.5 * np.add(thetapert_rough[j, :, :],
                                                  thetapert_rough[j - 1, :, :])
        wvelpert = wvels_list[i]

        slice_lev = np.where(np.abs(height - height_level) < 26)[0][0]

        wvelthetapert = np.multiply(wvelpert, thetapert)
        wvelperts_list.append(wvelpert[slice_lev, :, :])
        thetaperts_list.append(thetapert[slice_lev, :, :])

        wvelthetaperts_list.append(wvelthetapert)
        full_profs_list.append((wvelpert, thetapert, wvelthetapert))

    #flatten the arrays, TODO: make a function or class method
    wvelperts = np.array(wvelperts_list)
    thetaperts = np.array(thetaperts_list)
    [enum, ynum, xnum] = wvelperts.shape

    #wvelperts_slice = wvelperts[0]
    #thetaperts_slice = thetaperts[0]

    wvelperts = np.reshape(wvelperts, enum * ynum * xnum)
    thetaperts = np.reshape(thetaperts, enum * ynum * xnum)

    return wvelperts, thetaperts, full_profs_list, height, filenames


def write_h5(case_dict,h5_outfile='test.h5'):
    hvals_names = ['zg0', 'h', 'eltop_dthetadz', 'zf0', 'elbot_flux', 'h_flux',
                   'eltop_flux', 'deltatheta', 'mltheta', 'z1_GMa']
    hvals_columns = OrderedDict(zip(hvals_names, list(range(len(
        hvals_names)))))
    hvals_tup = make_tuple(hvals_columns)
    #scale order
    #defined in get_limits.py
    scale_names = ['rino', 'invrino', 'wstar', 'S', 'tau', 'mltheta',
                   'deltatheta', 'pi3', 'pi4', 'thetastar', 'c_delta']
    scale_columns = OrderedDict(zip(scale_names, list(range(len(
        hvals_names)))))
    scale_tup = make_tuple(scale_columns)


    avproflims_prefix = "/tera/users/nchaparr/"
    avproflims_dir = "/data/AvProfLims"
    invrinos_dir = "/data/invrinos"

    with h5py.File(h5_outfile, 'w') as f:
        the_dict = {}
        for casename,the_dict in case_dict.items():
            time_index = the_dict['time_index']
            start,step,stop = the_dict['time_list']
            dump_time_list, Times = Make_Timelists(start,step,stop)
            print('here are times: ', Times)
            #get heights and convective scales from text files
            hval_file = "{}{}{}".format(avproflims_prefix, casename,
                                        avproflims_dir)
            hvals = np.genfromtxt(hval_file)
            print('hvals: ', hvals)
            invrino_file = "{}{}{}".format(avproflims_prefix, casename,
                                           invrinos_dir)
            scales = np.genfromtxt(invrino_file)
            print('scales: ', scales)
            try:
                thetastar, wstar = scales[time_index, scale_tup.thetastar], scales[
                    time_index, scale_tup.wstar]
            except IndexError:
                print('skipping case {}: requested time index {}, total run is length {}'.format(
                    casename,time_index,scales.shape[0]))
                continue
            wvelperts, thetaperts, full_profs_list, height, sam_filelist = get_ensemble(
                casename, dump_time_list[time_index],
                hvals[time_index, hvals_tup.zg0])
            print(full_profs_list[0][0].shape)
            case = f.create_group(casename)
            keys = ['wvelpert', 'thetapert', 'wvelthetapert']
            for count, three_perts in enumerate(full_profs_list):
                run = case.create_group("{}".format(count))
                key_pairs = zip(keys, three_perts)
                for key, array in key_pairs:
                    print('writing ', case, key)
                    dset = run.create_dataset(
                        key, array.shape, dtype=array.dtype)
                    dset[...] = array[...]
            dset = case.create_dataset('height',
                                       height.shape,
                                       dtype=height.dtype)
            dset[...] = height[...]
            case.attrs['sam_filelist'] = json.dumps(sam_filelist)
            dset = case.create_dataset('scales',
                                       scales.shape,
                                       dtype=scales.dtype)
            dset[...] = scales[...]
            dset.attrs['scale_columns'] = json.dumps(scale_columns)
            dset.attrs['scale_file'] = invrino_file
            dset = case.create_dataset('hvals', hvals.shape, dtype=hvals.dtype)
            dset[...] = hvals[...]
            dset.attrs['height_columns'] = json.dumps(hvals_columns)
            dset.attrs['hval_file'] = hval_file
            dset = case.create_dataset('time_list',
                                       Times.shape,
                                       dtype=Times.dtype)
            dset.attrs['time_index'] = time_index
            dset[...] = Times[...]
        f['/'].attrs['case_dict']=json.dumps(case_dict,indent=4)
    return None


if __name__ == "__main__":

    import argparse, textwrap
    linebreaks = argparse.RawTextHelpFormatter
    descrip = textwrap.dedent(globals()['__doc__'])
    parser = argparse.ArgumentParser(formatter_class=linebreaks,
                                     description=descrip)
    parser.add_argument('-j', '--jfile', help='json file with run info', required=True)
    parser.add_argument('-r',
                        '--root',
                        help='root name of output hdf5',
                        required=True)
    args = parser.parse_args()


    with open(args.jfile,'r') as f:
        case_dict = json.load(f)
        outfile = 'profiles_{}.h5'.format(args.root)
        write_h5(case_dict,h5_outfile=outfile)


