# -*- coding: utf-8 -*-

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

from collections import Counter
#from neuropype_graph.utils_stats import compute_oneway_anova_fwe,compute_pairwise_ttest_fdr

def compute_signif_permuts(permut_df, permut_col = "Seed",session_col = "Session", start_col = 3, stop_col = 0):
    """
    args:
    compute significance of permutation over a df generated by gather_permuts 
    permut_df: original permutation results (pandas Dataframe)
    stop_col: last column to be included (in fast, excluded except if value is 0, in this case goes to the last column of the df
    
    return:
    all_p_higher,    all_p_lower: "vector of p_values obtained for 1 tail t-test in both direction, first session - second session"
    
    """
    
    print permut_df
    
    
    ################################################## check if permut_col exists and is consistent with permutation indexexing:
    seed_index = np.unique(permut_df[permut_col].values)
    
    print seed_index
    
    ### should start with -1
    assert seed_index[0] == -1,"Error, permut_col {} should start with -1".format(permut_col)
    
    expected_permut_indexes = range(len(seed_index)-1)
    
    print expected_permut_indexes
    
    ### should start at 0 and have all values in between
    assert all(x in seed_index[1:] for x in expected_permut_indexes),"Error, permut indexes should be consecutive and start with 0:  {} ".format(expected_permut_indexes)
    
    nb_permuts = len(expected_permut_indexes)
    
    print nb_permuts
    
    ### all unique values should have 2 different samples
    count_elements = Counter(permut_df[permut_col].values)
    
    print count_elements
    
    assert all(val == 2 for val in count_elements.values()), "Error, all permut indexes should have 2 and only 2 lines: {}".format(count_elements)
    
    ################################################## computing diff df
    
    #### orig diff
    #if stop_col == 0:
        #orig_df = permut_df[permut_df[permut_col] == -1].iloc[:,start_col:]
        
    #else:
        #orig_df = permut_df[permut_df[permut_col] == -1].iloc[:,start_col:stop_col]
        
    #print orig_df
    
    #orig_diff = orig_df.loc[0,].values - orig_df.loc[1,].values
    
    #print orig_diff
    
    #sign_orig_df = np.sign(orig_diff)
    
    #print sign_orig_df
    
    ### diff_df 
    
    if stop_col == 0:
        data_cols = permut_df.columns[start_col:]
        
    else:
        data_cols = permut_df.columns[start_col:stop_col]
        
    all_p_higher = np.zeros(shape = (len(data_cols)), dtype = 'float64') -1
    all_p_lower = np.zeros(shape = (len(data_cols)), dtype = 'float64') -1
                             
    for index_col,col in enumerate(data_cols):
        
        print col
        
        df_col = permut_df.pivot(index = permut_col, columns = session_col, values = col)
        
        print df_col
        
        df_col["Diff"] = df_col.iloc[:,0] - df_col.iloc[:,1]
        
        print df_col["Diff"]
        print df_col.shape
        
        if df_col["Diff"].iloc[0] > 0:
            sum_higher = np.sum((df_col["Diff"].iloc[1:] > df_col["Diff"].iloc[0]).values.astype(int))
            print "sum_higher:",sum_higher
            all_p_higher[index_col] = (sum_higher+1)/float(df_col.shape[0])
            
        else :
            sum_lower = np.sum((df_col["Diff"].iloc[0] > df_col["Diff"].iloc[1:]).values.astype(int))
            print "sum_lower:",sum_lower
            all_p_lower[index_col] = (sum_lower+1)/float(df_col.shape[0])
        
        #print df_col["Diff"] < df_col["Diff"][0]
    print all_p_higher
    print all_p_lower
    
    return all_p_higher,all_p_lower

def compute_mean_cormats(all_cormats,all_descriptors,descript_columns):

    print "In compute_mean_cormats"
    
    for column in descript_columns:
    
        assert column in all_descriptors.columns, "Error, {} not in {}".format(column,all_descriptors.columns)

    dict_mean = {}
    
    for elem, lines in all_descriptors.groupby(by = descript_columns):
    
        print elem
        #print lines
        print lines.index
        
        print all_cormats.shape
        
        elem_cormats = all_cormats[lines.index,:,:]
        
        print elem_cormats.shape
        
        mean_elem = np.mean(elem_cormats,axis = 0)
        
        print mean_elem.shape
        
        dict_mean[elem] = mean_elem
        
    return dict_mean




def compute_stats_cormats(all_cormats,all_descriptors,descript_columns, groups = []):

    print all_cormats.shape
    
    for column in descript_columns:
    
        assert column in all_descriptors.columns, "Error, {} not in {}".format(column,all_descriptors.columns)

    dict_stats = {}
    
    for column in descript_columns:
        
        if len(groups) == 0:
            groups = all_descriptors[column].unique().tolist()
        
        ############## compute F-test over matrices
        
        list_of_list_matrices = [all_cormats[all_descriptors[all_descriptors[column] == cond_name].index,:,:] for cond_name in groups]
        
        print list_of_list_matrices
        print len(list_of_list_matrices)
        print np.array(list_of_list_matrices).shape
        
        
        signif_adj_mat = compute_oneway_anova_fwe(list_of_list_matrices,cor_alpha = 0.05, uncor_alpha = 0.01)
        
        print signif_adj_mat
             
        dict_stats["F-test"] = signif_adj_mat
        
        for combi_pair in combinations(groups,2):
            pair_name = "-".join(combi_pair)
            print pair_name
            
            print combi_pair[0]
            
            try:
                signif_signif_adj_mat = compute_pairwise_ttest_fdr(X = list_of_list_matrices[groups.index(combi_pair[0])],
                                                               Y = list_of_list_matrices[groups.index(combi_pair[1])],
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
	