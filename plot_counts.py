from matplotlib import pyplot as plt
import numpy as np
import sys
import pandas as pd
ndtime=sys.argv[1]
infile='all_{}.h5'.format(ndtime)
with pd.HDFStore(infile) as store:
    df_all=store['/all_cases']

print(df_all.columns)

cases = set(df_all['case'])
columns=df_all.columns
quad_counts=[item for item in columns if 'Wcount' in item]
print(quad_counts)
plt.close('all')
for count,case in enumerate(cases):
    fig, ax = plt.subplots(1,1)
    print('plotting ',count,case)
    df_case=df_all.loc[df_all['case'] == case]
    height=df_case['zenc_nd_height'].values
    the_sum=np.zeros_like(height)
    for curve in quad_counts:
        counts=df_case[curve].values
        the_sum=the_sum+counts
        ax.plot(counts,height,label=curve)
    ax.legend(loc='best')
    ax.set(title='{}, ndtime: {} -- quadrant counts'.format(case,ndtime))
    filename='quad_{}.png'.format(ndtime)
    fig.savefig(filename)
plt.show()    
    
