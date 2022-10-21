#------------------------------------------------------------------
# Exhaustive TURF Analysis
# author: e.grant
#------------------------------------------------------------------

# 1.Import needed libraries
import pandas as pd
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import ttest_rel
import os as os
cwd = os.getcwd()
if not os.path.exists('output'):
    os.makedirs('output')

from datetime import datetime 
startTime = datetime.now()
 
###################################################################################################
###################################################################################################

   # A. Define Parameters
proj = 'FXXXX'    #------------------------- Project Number
filename = 'utilities.csv' #------------------ Name of the utilities file
n_kv = 44    #-------------------------------- Total number of key visuals / concepts    
k = 2       # ----------------------------------- Number of elements to be drawn
n_fixed = 0  # ------------------------------- How many fixed values we will have (if none = 0)
#fixed1 = 9   # ----------------------------- Copy below as many times as needed and modify "Defining fixed values code"     
#fixed2 = 14   # ----------------------------- Copy below as many times as needed and modify "Defining fixed values code"     
#fixed3 = 18   # ----------------------------- Copy below as many times as needed and modify "Defining fixed values code"     
#fixed4 = 7
#fixed5 = 11
#fixed6 = 12
#fixed7 = 10
#fixed8 = 21
#fixed9 = 6
#fixed10 = 8
#fixed11 = 15
#fixed12 = 16
#fixed13 = 17
#fixed14 = 18
#fixed15 = 19
#fixed16 = 20
#fixed17 = 21

prod_wei = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
prod_wei2 = pd.DataFrame(list(prod_wei))
prod_wei2 = pd.DataFrame.transpose(prod_wei2)

sigtest = 0.05    # -------------------------- Significance either (0.05 or 0.1)
corrmap = 'n'     # -------------------------- Generate correlation map? options = y/n  
sourcing = False  # -------------------------- Use sourcing analysis? no = False

# B. Here we define the buckets (Exhaustive list every concept MUST be in one and only one bucket)
comp = [44] # --- Competitor bucket Concepts / Key Visuals
n_buckets = 2   # ---------------------------- number of client'sbuckets
buk_1 = [1, 2, 3, 4,5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21]  # ---------------------- Bucket name 1
buk_2 = [22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43] # --- Bucket name 2

# C. Capping of buckets one line per bucket (Maximun number of concepts per bucket allowed - default = k)
lim_buk1 = 21    # ---------------------------- Capping 
lim_buk2 =  22   # ---------------------------- Capping 

# D. Here we as many restrictions as we want (buckets non-allowed to combine)
n_rest = 0      # --------------------------- Total number of restrictions 
# Example: rest_1 = [buk_1,buk_2]

# E. Parameter for significant testing
sigtest = 0.05  # --------------------------- Significance either (0.05 or 0.1)

# F. Pre-work (Do not change!)
lst = list(range(1, n_kv+1))
head = list(range(1, k+1))
head.append('None')
cte = list(set(lst) - set(comp))
lst_cte = ['KV' + str(x) for x in cte]
cte[:] = [x - 1 for x in cte]
len_com = len(comp)
buk = []
lim = []
for i in range(1,n_buckets+1):
    buk.append(locals()['buk_' + str(i)])
    lim.append(locals()['lim_buk' + str(i)])

###################################################################################################
###################################################################################################

# 1. Read Utilities
util = pd.read_csv(filename) 

# 1a. Keep breaks as arrays and delete
wei = pd.DataFrame(util['wei'])

# Delete overhead from the dataframe!
del util['wei']
del util['id']

# Exponentiation of the utilities
util = np.exp(util)

# 2. Correlation matrix of all key visuals
if (corrmap == 'y'):
    corr = util.corr()
    corr.values[[np.arange(n_kv)]*2] = np.nan
    corr = corr.fillna('')
    corr.to_csv('output/CorrMatrix.csv')
    plt.figure(figsize=(20,8))
    sns_plot = sns.heatmap(util.corr(),cmap="PiYG",annot=True,vmin=-1, vmax=1)
    fig = sns_plot.get_figure()
    fig.savefig('output/'+ proj + '_CorrHeatMap', dpi=500)

# Moving on!!!
util = pd.read_csv(filename)
del util['id']
del util['wei']
util = np.exp(util)

# 3. Matrix of combinations
elem = np.arange(min(min(buk_1),min(buk_2)), max(max(buk_1),max(buk_2))+1)
nck = pd.DataFrame(list(combinations(elem,k)))
#nck_bin = pd.get_dummies(nck.stack()).sum(level=0)
nck_bin = pd.get_dummies(nck.stack()).groupby(level=0).sum()

#nck_bin[comp] = 1
for x in range(min(comp),max(comp)+1):
    nck_bin[x] = 1

# 4. Clean the Matrix to keep only relevant combinations
nck_bin['sum'] = nck_bin.sum(axis=1)    
nck_bin.drop(nck_bin[nck_bin['sum'] != (len_com+k)].index, inplace=True)

# 3.1 Uncomment and add as many lines as neeed
#if (n_fixed != 0):
#    nck_bin.drop(nck_bin[nck_bin[fixed1] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed2] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed3] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed4] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed5] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed6] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed7] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed8] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed9] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed10] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed11] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed12] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed13] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed14] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed15] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed16] != 1].index, inplace=True)
#    nck_bin.drop(nck_bin[nck_bin[fixed17] != 1].index, inplace=True)

# 4. Functions to limit combinations among and within buckets (each bucket can have diff cap)
def bucket(bucket,binary):
    const = pd.DataFrame(list(combinations(bucket,2)))
    const = pd.get_dummies(const.stack()).sum(level=0)
    const=pd.merge(nck_bin, const, on=bucket, how='inner')
    binary=pd.concat([binary, const]).drop_duplicates(keep=False)
    binary=pd.concat([binary, const])
    return (binary)

def constraints(bin_matrix,nbuk,bucket,limit):
    bin_matrix['Sum_B'+str(nbuk)] = bin_matrix[bucket].sum(axis=1)
    bin_matrix.drop(bin_matrix[bin_matrix['Sum_B'+str(nbuk)] > limit].index, inplace=True)
    del nck_bin['Sum_B'+str(nbuk)]
    return (bin_matrix)

# 5. Capping of combinations within buckets
for i in range(1,n_buckets+1):
    nck_bin = constraints(nck_bin,1,buk[i-1],lim[i-1])

# 6. Capping of combinations between buckets (non-allowed)
if (n_rest != 0):
    all_const = [(x,y) for x in buk_1 for y in buk_2]
    for i in range(1,len(all_const)+1):
        const = list(all_const[i-1])
        nck_bin['const'] = nck_bin[const].sum(axis=1)
        nck_bin.drop(nck_bin[nck_bin['const'] == 2].index, inplace=True)
        del nck_bin['const'] 

# 9. Loop over all combinations and save SoP results
n_comb = nck_bin.shape[0]
del nck_bin['sum']
nck_bin['None'] = 1

# Empty data frames for storing results
ave_sop2 = pd.DataFrame()
sum_SoP   = pd.DataFrame()
KVs       = pd.DataFrame()
Esc       = pd.DataFrame()

for i in range(1,n_comb+1):    
    pos = i
    nck_bin2 = pd.DataFrame(nck_bin.iloc[pos-1:pos,:])  
    kv = nck_bin2.apply(lambda row: [f"KV{row.name}" if x == 1 else x for x in row], axis=0)
    kv = kv.loc[:, (kv != 0).any(axis=0)]
    kv.drop(comp, axis=1, inplace=True)
    kv.columns = head
    del kv['None']
    KVs = pd.concat([KVs,kv])
    
    nck_bin3 = pd.DataFrame(nck_bin2.values * prod_wei2.values)
    util2 = pd.DataFrame(util.values * nck_bin3.values, columns=util.columns)
    util2['total'] = util2.sum(axis=1)
    exp = util2.div(util2.total, axis='index')
    exp['esc'+str(pos)] = exp[lst_cte].sum(axis=1)
    exp = pd.DataFrame(exp.values * wei.values, columns=exp.columns)    
    t_Esc = exp['esc'+str(pos)]
    Esc = pd.concat([Esc,t_Esc], axis = 1) 
    ave_sop = pd.DataFrame(exp.mean(), columns=[i])
    data = ave_sop.iloc[cte,:].sum()
    kvs = ave_sop.iloc[cte,:]
    data2 = ave_sop.iloc[:,:]
    sum_SoP = pd.concat([sum_SoP,data])
    ave_sop2 = pd.concat([ave_sop2,ave_sop], axis = 1,sort=False) 
    
# 11. Tidy up results for saving - No one likes messy results
KVs.index = np.arange(1,len(KVs)+1)
sum_SoP = pd.concat([sum_SoP,KVs], axis = 1)
sum_SoP = sum_SoP.sort_values(by=0, ascending = False)
sum_SoP=sum_SoP.rename(columns = {0:'SoP'})

# 11. Paired T Test for evaluation of turf performance
l_esc = sum_SoP.index.tolist()
p_value = [1]
ctrl = 0
sig = ["=="]
for e in range (1,len(l_esc)):
    sol = ttest_rel(Esc['esc'+str(l_esc[ctrl])], Esc['esc'+str(l_esc[e])])
    p_value.append(sol[1])
    if sol[1] < sigtest:
        ctrl = e
        sig.append("***")
    else:
        sig.append("==")

# 12. Let's save our neat results
sum_SoP["p_value"] = p_value
sum_SoP["sigtest"] = sig
sum_SoP.to_csv('output/'+ proj + '_TURF'+str(k)+'_Results.csv')
Esc.to_csv('output/'+ proj + '_TURF'+str(k)+'_SOP.csv')
ave_sop2.to_csv('output/'+ proj + '_TURF'+str(k)+'_Avg_SOP.csv')
#exp.to_csv('output/'+ proj + '_TURF' +str(k)+'_ForCalibration.csv')

print("------------------------------------------------")
print("Tool evaluated " +str(n_comb) + " combinations")
print("Process for " +str(k) + " concepts took: " + str(datetime.now() - startTime))
print("------------------------------------------------")
 