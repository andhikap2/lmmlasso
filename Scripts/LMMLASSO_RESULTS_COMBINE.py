import pandas as pd
import glob

path = r'/srv/uom-data1-q.unimelb.edu.au/6300-afournier/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Summary_Results' # use your path
all_files = glob.glob(path + "/*[0-9]*.csv*")


li = []

for filename in yeah:
    df = pd.read_csv(filename,index_col=None)
    a=df.iloc[0,:]
    li.append(a)
   

    #li=np.concatenate((li,a),axis=1)
    li=pd.concat((li,a),axis=0)

frame = pd.concat(li, axis=1, ignore_index=True)
frame=frame.transpose()
frame=frame.iloc[:,1:9]
frame.to_csv('LMMLASSO_RESULTS_COMBINED.csv',frame)