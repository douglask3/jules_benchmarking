import os
from os.path import isfile 
from   io     import StringIO
import numpy  as np
import pandas as pd
import csv

import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy import stats

from pdb import set_trace as browser

file = "data/JULES-ES-UKESM1 - Tidied.csv"

def binnize(x, y, xmin, xmax = None, xbin_size = 1):
    if xmax is None:
        xmax = xmin[1] if len(xmin) == 2 else xmin[0] + xbin_size
        xmin = xmin[0]

    return y[np.logical_and(x >= xmin, x <= xmax)]
    
def bin_and_ttest_series(x, y, hist_cut, xbin_size, xwind_size, z = None):
    x = np.array(x).astype(float)
    y = np.array(y).astype(float)

    hist_y = binnize(x, y, hist_cut)
    
    def fun(xs): return binnize(x, y, [xs], xbin_size = xbin_size)

    index = np.arange(min(x), max(x)- xbin_size, xwind_size)
    comp_y = [fun(xs) for xs in index]

    def fun(cy): return stats.ttest_ind(hist_y, cy, alternative = 'greater').pvalue
    pvs = np.array([fun(cy) for cy in comp_y])
    
    for i in range(len(pvs)): pvs[i] = max(pvs[i:])     
    

    sigPnt = np.where(pvs < 0.01)[0][0]
    sigX = x[sigPnt]+ xbin_size
    
    if z is not None: sigZ = [zi[np.argmin(abs(x - sigX))] for zi in z]
    return pvs, sigX, sigZ


dat = pd.read_csv(file)


### example
time = dat["Variable"][2:]
cover = dat["Tree Cover"][2:]

tas = dat['Tas'][2:]

pvs, sigYr, sigTas = bin_and_ttest_series(time, cover, [1995, 2015], 20, 1, [tas])

regions =  set(dat.iloc[0,1:])
ssps = set(dat.iloc[1,1:])

def forRegion(region, ssp):
    print(region)
    print(ssp)
    if isinstance(ssp, str) and region != 'global': 
    
        def find_row(ssp, region, var):
            varID = np.array([var in i for i in np.array(dat.keys())])
            id = np.logical_and(dat.iloc[1,:] == ssp, dat.iloc[0,:] == region)
            id = np.logical_and(id, varID)
           
            try:
                id = np.where(id)[0][0]    
            except:
                browser()        
            return dat.iloc[2:, id]
        
        tasrg = find_row(ssp, region, 'Tas')
        
        cover = find_row(ssp, region, "Tree Cover")        
        tasgl = find_row(ssp, "global", 'Tas')
        
        #id = np.logical_and(dat.iloc[1,:] == ssp, dat.iloc[0,:] == "global")
        #id = np.where(id)[0][0]
        #tas = dat.iloc[2:, id]
        
        pvs, sigYr, sigTas = bin_and_ttest_series(time, cover, [1995, 2016], 20, 1, 
                                                 [tasgl, tasrg])
        if region == "Aus": browser()
        sigTas = np.round(np.array(sigTas, dtype = 'float'), 2)
        out = str(sigYr)[0:4] 
        for i in sigTas: out = out + ', ' + str(i)
        
        return out
    else: 
        return '--'

out = np.array([[forRegion(region, ssp) for region in regions] for ssp in ssps])
out = pd.DataFrame(out.T, columns = ssps, index = regions)
out.to_csv('tree_divergance_time_GW.csv')
browser()


