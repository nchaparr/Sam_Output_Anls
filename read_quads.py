import pandas as pd
import glob,re
import numpy as np

date_list = ["Mar52014", "Jan152014_1", "Dec142013", "Nov302013", 
             "Mar12014", "Dec202013", "Dec252013"] #"","", 

vars=['upwarm_bar','downwarm_bar','upcold_bar', 
         'downcold_bar','wvelthetapert_bar']

heightfile='/newtera/tera/phil/nchaparr/python/Plotting/Mar12014/data/heights0000000600'
heights=np.genfromtxt(heightfile)

columns=['upwarm_bar','downwarm_bar','upcold_bar', 
         'downcold_bar','wvelthetapert_bar']
timeval=re.compile(".*format(\d{10,10})")

df_dict={}
with pd.HDFStore('quads.h5','w') as store:
    for case in date_list:
        rootdir="/home/phil/repos/Sam_Output_Anls/dump/{}/data".format(case)
        files=glob.glob("{}/flux_quads*".format(rootdir))
        for the_file in files:
            the_match=timeval.match(the_file)
            num=the_match.groups(1)[0]
            the_time=int(num)
            numbers=np.genfromtxt(the_file)
            df=pd.DataFrame(numbers,columns=columns)
            nodename='/{}/t_{}'.format(case,the_time)
            store.put(nodename,df,format='table')
    store.put('/heights',pd.Series(heights))

