"""
given paper_table.h5 produced by gm_numbers.py containing the
run characteristics and a target non-dimensional time, write
a json file with the time idices of each case that produce a timestep closest
to that non-dimnesional time

example:  python calc_index.py -i paper_table.h5 -t 22.8

"""

import json
import pandas as pd
from Make_Timelist import Make_Timelists
import numpy as np

def get_dict(case_file,target_time):
    
    with pd.HDFStore(case_file,'r') as store:
        df_table = store['cases']

    run_dict={}    
    for index, row in df_table.iterrows():
        line_dict = row.to_dict()
        case = line_dict['name']
        start, step, stop = 1, 600, 28800
        if case == 'Nov302013':
            step = 900
        dump_time_list, Times = Make_Timelists(start, step, stop)
        nd_times = (Times*3600.)*line_dict['N']
        #
        # calc zenc/L0
        #
        zenc_l0_times=np.sqrt(2*nd_times)
        print(nd_times, line_dict['N'])
        time_index=int(np.searchsorted(zenc_l0_times,time_nd_target))
        try:
            print('case: {}, target: {}, closest ndtime: {}'.format(case,time_nd_target,nd_times[time_index]))
        except IndexError:
            print('Dropping case: {}, target: {}, last nd_time: {}, L0: {}'.format(case,time_nd_target,nd_times[-1],line_dict['L0']))
            continue
        run_dict[case] = line_dict
        #print(run_dict[case]['legends'])
        run_dict[case]['time_list'] = (start,step,stop)
        time_index = int(time_index)
        run_dict[case]['time_index']=time_index
        print(time_index, case, len(nd_times))
        run_dict[case]['nd_time'] = nd_times[time_index]
        run_dict[case]['zenc_l0'] = zenc_l0_times[time_index]
    return run_dict
         

if __name__ == "__main__":

    import argparse, textwrap
    linebreaks = argparse.RawTextHelpFormatter
    descrip = textwrap.dedent(globals()['__doc__'])
    parser = argparse.ArgumentParser(formatter_class=linebreaks,
                                     description=descrip)
    parser.add_argument('-i', '--input', help='pandas table file with run information', required=True)
    parser.add_argument('-t', '--nd_time', help='non-dimensional target time',type=float, required=True)
    args = parser.parse_args()

    
    infile = args.input
    time_nd_target = args.nd_time
    run_dict = get_dict(infile,time_nd_target)

    outfile = 'new_index_list_time_nd_{}.json'.format(time_nd_target)
    with open(outfile,'w') as f:
        print(run_dict)
        json.dump(run_dict,f,indent=4)

    
