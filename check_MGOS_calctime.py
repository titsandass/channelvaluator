import pickle, os

MGOScalctimePath = os.getcwd() + "\\channels\\MGOS_calctime_NOQTF.pickle"
MGOScalctimePath2 = os.getcwd() + "\\channels\\MGOS_calctime.pickle"

with open(MGOScalctimePath, 'rb') as f:
    calcTime = pickle.load(f)

# with open(MGOScalctimePath2, 'rb') as f:
#     calcTime2 = pickle.load(f)

pass