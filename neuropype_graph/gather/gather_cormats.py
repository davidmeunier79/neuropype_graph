
import os

import pandas as pd
import numpy as np
#import statsmodels as stat
#import statsmodels.formula.api as smf
#import statsmodels.api as sm
#import matplotlib.pyplot as plt
#import nibabel.gifti as gio
#from statsmodels.stats.outliers_influence import OLSInfluence

from itertools import product,combinations
from numpy import isnan, nan, logical_not, logical_or

def isInAlphabeticalOrder(word):
    return list(word) == sorted(word)
    #return word==''.join(sorted(word))

def return_all_iter_cormats(cormat_path ,iterables ,iternames ):
       
    print zip(*iterables)
    
    all_iter_cormats = []
    all_descriptors = []
    
    assert isInAlphabeticalOrder(iternames), "Warning, iternames are not in alphabetical oroder, check the iterables order as well"
    
    for iter_obj in product(*iterables):
        
        print iter_obj
        
        assert len(iter_obj) == len(iternames), "Error, different number of iternames and iterables"
        
        iter_dir = "".join(["_" + zip_iter[0] + "_" + zip_iter[1] for zip_iter in zip(iternames,iter_obj)])
                            
        print iter_dir
        
        cormat_file = os.path.join(cormat_path,iter_dir,"compute_conf_cor_mat","Z_cor_mat_resid_ts.npy")
        
        if os.path.exists(cormat_file):
            cormat = np.load(cormat_file)
            print cormat.shape
            
            all_iter_cormats.append(cormat)
            all_descriptors.append(iter_obj)
            
        else:
            print "Warning, file {}  could not be found".format(cormat_file)
        
    return np.array(all_iter_cormats),pd.DataFrame(all_descriptors,columns = iternames)
    

if __name__ =='__main__':
	
	test1 = isInAlphabeticalOrder(["a","b","c"])
	print test1
	
	test2 = isInAlphabeticalOrder(["ab","ba","ca"])
	print test2
	