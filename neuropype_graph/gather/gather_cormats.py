
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

from neuropype_graph.utils_stats import compute_oneway_anova_fwe,compute_pairwise_ttest_fdr

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

def compute_mean_cormats(all_cormats,all_descriptors,descript_columns):

    for column in descript_columns:
    
        assert column in all_descriptors.columns, "Error, {} not in {}".format(column,all_descriptors.columns)

    dict_mean = {}
    
    for elem, lines in all_descriptors.groupby(by = descript_columns):
    
        print elem
        print lines
        print lines.index
        
        elem_cormats = all_cormats[lines.index]
        
        print elem_cormats.shape
        
        mean_elem = np.mean(elem_cormats,axis = 0)
        
        print mean_elem.shape
        
        dict_mean[elem] = mean_elem
        
    return dict_mean

def compute_stats_cormats(all_cormats,all_descriptors,descript_columns):

    print all_cormats.shape
    
    for column in descript_columns:
    
        assert column in all_descriptors.columns, "Error, {} not in {}".format(column,all_descriptors.columns)

    dict_stats = {}
    
    for column in descript_columns:
        
        cond_names = all_descriptors[column].unique().tolist()
        
        ############## compute F-test over matrices
        
        list_of_list_matrices = [all_cormats[all_descriptors[all_descriptors[column] == cond_name].index,:,:] for cond_name in cond_names]
        
        print list_of_list_matrices
        print len(list_of_list_matrices)
        print np.array(list_of_list_matrices).shape
        
        
        signif_adj_mat = compute_oneway_anova_fwe(list_of_list_matrices,cor_alpha = 0.05, uncor_alpha = 0.01)
        
        print signif_adj_mat
             
        dict_stats["F-test"] = signif_adj_mat
        
        for combi_pair in combinations(cond_names,2):
            pair_name = "-".join(combi_pair)
            print pair_name
            
            print combi_pair[0]
            
            try:
                signif_signif_adj_mat = compute_pairwise_ttest_fdr(X = list_of_list_matrices[cond_names.index(combi_pair[0])],
                                                               Y = list_of_list_matrices[cond_names.index(combi_pair[1])],
                                                               cor_alpha = 0.05, uncor_alpha = 0.01,paired = True,old_order = False)
                
                print signif_signif_adj_mat
                dict_stats["T-test_" + pair_name] = signif_signif_adj_mat
                
            except AssertionError:
                print "Stop running after {} was wrong".format(pair_name)
                
    return dict_stats

if __name__ =='__main__':
	
	test1 = isInAlphabeticalOrder(["a","b","c"])
	print test1
	
	test2 = isInAlphabeticalOrder(["ab","ba","ca"])
	print test2
	