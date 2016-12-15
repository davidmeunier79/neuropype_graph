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

def compute_rada_df(iter_path,df):

    from neuropype_graph.utils_mod import get_modularity_value_from_lol_file
    from neuropype_graph.utils_mod import get_values_from_global_info_file
    from neuropype_graph.utils_mod import get_path_length_from_info_dists_file
        
    ########### modularity
    
    modularity_file = os.path.join(iter_path,"community_rada","Z_List.lol")

    print modularity_file
    

    
    
    if os.path.exists(modularity_file):
    
        mod_val = get_modularity_value_from_lol_file(modularity_file)
    #else:
        #mod_val = np.nan
            
        print mod_val
        
        df['Modularity'] = mod_val
        
    print df
    
    #################### info_global 
    
    global_info_file = os.path.join(iter_path,"net_prop","Z_List-info_global.txt")
    
    print global_info_file
    
    if os.path.exists(global_info_file):
    
    
        global_info_values = get_values_from_global_info_file(global_info_file)
        
        print global_info_values
        
        df.update(global_info_values)
        
        print df
        
    ##################### info_dists
    
    path_length_file = os.path.join(iter_path,"net_prop","Z_List-info_dists.txt")

    print path_length_file
    
    if os.path.exists(path_length_file):
    
        mean_path_length,diameter,global_efficiency = get_path_length_from_info_dists_file(path_length_file)
        
        print mean_path_length,diameter
        
        df['Mean_path_length'] = str(mean_path_length)
        df['Diameter'] = str(diameter)
        df['Global_efficiency'] = str(global_efficiency)
    
    print df
            

def compute_nodes_rada_df(iter_path,df):

    from neuropype_graph.utils_mod import get_modularity_value_from_lol_file
    from neuropype_graph.utils_mod import get_values_from_global_info_file
    from neuropype_graph.utils_mod import get_path_length_from_info_dists_file
        
    ########### modularity
    
    modularity_file = os.path.join(iter_path,"community_rada","Z_List.lol")

    print modularity_file
    

    
    
    if os.path.exists(modularity_file):
    
        mod_val = get_modularity_value_from_lol_file(modularity_file)
    #else:
        #mod_val = np.nan
            
        print mod_val
        
        df['Modularity'] = mod_val
        
    print df
    
    #################### info_global 
    
    global_info_file = os.path.join(iter_path,"net_prop","Z_List-info_global.txt")
    
    print global_info_file
    
    if os.path.exists(global_info_file):
    
    
        global_info_values = get_values_from_global_info_file(global_info_file)
        
        print global_info_values
        
        df.update(global_info_values)
        
        print df
        
    ##################### info_dists
    
    path_length_file = os.path.join(iter_path,"net_prop","Z_List-info_dists.txt")

    print path_length_file
    
    if os.path.exists(path_length_file):
    
        mean_path_length,diameter,global_efficiency = get_path_length_from_info_dists_file(path_length_file)
        
        print mean_path_length,diameter
        
        df['Mean_path_length'] = str(mean_path_length)
        df['Diameter'] = str(diameter)
        df['Global_efficiency'] = str(global_efficiency)
    
    print df
            

def compute_nodes_rada_df(local_dir,list_df,gm_coords,coords_file,labels_file)

    from neuropype_graph.utils_net import read_lol_file,read_Pajek_corres_nodes
    from neuropype_graph.utils_dtype_coord import where_in_coords
    
    #### Z_List
    Pajek_file = os.path.join(local_dir,"prep_rada","Z_List.net")
    
    if os.path.exists(coords_file) and os.path.exists(Pajek_file) and os.path.exists(labels_file):
                        
        #### labels
        labels = np.array([line.strip() for line in open(labels_file)], dtype = str)
        
        #### MNI coordinates
        coords = np.array(np.loadtxt(coords_file),dtype = int)
        
        print coords.shape
            
        #### nodes in the connected graph
        
        node_corres = read_Pajek_corres_nodes(Pajek_file)
        
        print np.min(node_corres),np.max(node_corres)
        
        print node_corres.shape
            
                
        
        ### node_coords
        node_coords = coords[node_corres,:]
        
        print node_coords.shape
        
        node_labels = labels[node_corres,:].reshape(-1,1)
        print node_labels.shape
        
        
        ### where_in_gm_mask 
        where_in_gm_mask = where_in_coords(node_coords,gm_coords)
        
        where_in_gm_mask = where_in_gm_mask.reshape(where_in_gm_mask.shape[0],1)
        
        #print where_in_gm_mask
        print where_in_gm_mask.shape
        
        list_df.append(pd.DataFrame(np.concatenate((where_in_gm_mask,node_labels,node_coords),axis = 1),columns = ['Where_in_GM_mask','labels','MNI_x','MNI_y','MNI_z']))
            
    #### info nodes
    info_nodes_file = os.path.join(local_dir,"net_prop","Z_List-info_nodes.txt")
    
    if os.path.exists(info_nodes_file) :
        
        ## loading info_nodes
        df_node_info = pd.read_table(info_nodes_file)
        
        print "Info nodes:" 
        print df_node_info.shape
        
        list_df.append(df_node_info)
        
    else:
        print "Info nodes not found:" 
        print info_nodes_file
    
    #### modules /community_vect    
    partition_file = os.path.join(local_dir,"community_rada","Z_List.lol")
    
    if os.path.exists(partition_file) :
        
        
        ##loading partition_file
        community_vect = read_lol_file(partition_file)
        
        print "community_vect:" 
        print community_vect.shape
            
        list_df.append(pd.DataFrame(community_vect,columns = ['Module']))
        
    #### node roles    
    roles_file = os.path.join(local_dir,"node_roles","node_roles.txt")
    
    part_coeff_file = os.path.join(local_dir,"node_roles","all_participation_coeff.txt")
    
    Z_com_degree_file = os.path.join(local_dir,"node_roles","all_Z_com_degree.txt")
    
    if os.path.exists(roles_file) and os.path.exists(part_coeff_file) and os.path.exists(Z_com_degree_file):
        
        #### loding node roles
        
        node_roles = np.array(np.loadtxt(roles_file),dtype = int)
        
        print "node_roles:" 
        print node_roles.shape
        
        
        print "part_coeff:" 
        part_coeff = np.loadtxt(part_coeff_file)
        
        part_coeff = part_coeff.reshape(part_coeff.shape[0],1)
        
        print part_coeff.shape
        
        
        
        print "Z_com_degree:" 
        Z_com_degree = np.loadtxt(Z_com_degree_file)
        
        Z_com_degree = Z_com_degree.reshape(Z_com_degree.shape[0],1)
        
        print Z_com_degree.shape
        
        list_df.append(pd.DataFrame(np.concatenate((node_roles,part_coeff,Z_com_degree),axis = 1),columns = ['Role_quality','Role_quantity','Participation_coefficient','Z_community_degree']))
        

def compute_signif_permuts(permut_df, permut_col = "Seed",session_col = "Session", start_col = 0, stop_col = 0):
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
    
    if stop_col == 0:
        data_cols = permut_df.columns[start_col:]
        
    else:
        data_cols = permut_df.columns[start_col:stop_col]
        
    print data_cols
    
    all_p_higher = np.zeros(shape = (len(data_cols)), dtype = 'float64') -1
    all_p_lower = np.zeros(shape = (len(data_cols)), dtype = 'float64') -1
    
    cols = []
    
    for index_col,col in enumerate(data_cols):
        
        
        print index_col,col
        
        #print permut_df
        
        print permut_col,session_col
        
        df_col = permut_df.pivot(index = permut_col, columns = session_col, values = col)
        
        print df_col
        
        
        df_col["Diff"] = pd.to_numeric(df_col.iloc[:,0]) - pd.to_numeric(df_col.iloc[:,1])
        
        print df_col["Diff"]
        print df_col.shape
        
        if df_col["Diff"].iloc[0] > 0:
            sum_higher = np.sum((df_col["Diff"].iloc[1:] > df_col["Diff"].iloc[0]).values.astype(int))
            print "sum_higher:",sum_higher
            all_p_higher[index_col] = (sum_higher+1)/float(df_col.shape[0])
            
        elif df_col["Diff"].iloc[0] < 0 :
            sum_lower = np.sum((df_col["Diff"].iloc[0] > df_col["Diff"].iloc[1:]).values.astype(int))
            print "sum_lower:",sum_lower
            all_p_lower[index_col] = (sum_lower+1)/float(df_col.shape[0])
        
        else :
            print "not able to do diff"
            
        cols.append(col)
        
        #print df_col["Diff"] < df_col["Diff"][0]
    print all_p_higher
    print all_p_lower
    
    print all_p_higher, all_p_lower
        
    df_res = pd.DataFrame([all_p_higher, all_p_lower],columns= cols)
    df_res.index = ["Higher","Lower"]
        
    print df_res
    
    return df_res

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
	