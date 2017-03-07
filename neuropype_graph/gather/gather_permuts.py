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
            

def compute_nodes_rada_df(local_dir,gm_coords,coords_file,labels_file):

    from neuropype_graph.utils_net import read_lol_file,read_Pajek_corres_nodes
    from neuropype_graph.utils_dtype_coord import where_in_coords
    
    import glob
    
    #### Z_List
    Pajek_files = glob.glob(os.path.join(local_dir,"prep_rada","*.net"))
    
    assert len(Pajek_files) == 1, "Error, no .net file found in {} prep_rada".format(local_dir)
    
    if len(Pajek_files) == 1:
        
        Pajek_file = Pajek_files[0]
        
    list_df = []
            
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
        
        node_labels = labels[node_corres].reshape(-1,1)
        print node_labels.shape
        
        
        ### where_in_gm_mask 
        where_in_gm_mask = where_in_coords(node_coords,gm_coords)
        
        where_in_gm_mask = where_in_gm_mask.reshape(where_in_gm_mask.shape[0],1)
        
        #print where_in_gm_mask
        print where_in_gm_mask.shape
        
        list_df.append(pd.DataFrame(np.concatenate((where_in_gm_mask,node_labels,node_coords),axis = 1),columns = ['Where_in_GM_mask','labels','MNI_x','MNI_y','MNI_z']))
    else:
        print "Missing {},{} or {}".format(Pajek_file,coords_file,labels_file)
        
    #### info nodes
    info_nodes_files = glob.glob(os.path.join(local_dir,"net_prop","*-info_nodes.txt"))
    
    if len(info_nodes_files) == 1:
        
        info_nodes_file = info_nodes_files[0]
     
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
    partition_files = glob.glob(os.path.join(local_dir,"community_rada","*.lol"))
    
    if len(partition_files) == 1:
        
        partition_file = partition_files[0]
     
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
	
    return list_df
######################################## computing permutation-based stats per nodes, over several sheetnames #########################################
 
def compute_signif_permuts(permut_df, permut_col = "Seed",session_col = "Session", start_col = 0, stop_col = 0, columns = []):
    """
    args:
    compute significance of permutation over a df generated by gather_permuts 
    permut_df: original permutation results (pandas Dataframe)
    stop_col: last column to be included (in fact, excluded except if value is 0, in this case goes to the last column of the df
    
    return:
    all_p_higher,    all_p_lower: "vector of p_values obtained for 1 tail t-test in both direction, first session - second session"
    
    """
    
    #print permut_df
    
    
    ################################################## check if permut_col exists and is consistent with permutation indexexing:
    seed_index = np.unique(permut_df[permut_col].values)
    
    print seed_index
    
    ### should start with -1
    if not seed_index[0] == -1:
        
        print "Error, permut_col {} should start with -1".format(permut_col)
        #0/0
        return pd.DataFrame()
    
    expected_permut_indexes = range(len(seed_index)-1)
    
    print expected_permut_indexes
    
    ### should start at 0 and have all values in between
    if not all(x in seed_index[1:] for x in expected_permut_indexes):
        print "Error, permut indexes should be consecutive and start with 0:  {} ".format(expected_permut_indexes)
        
        #0/0
        
        #return pd.DataFrame()
    
    nb_permuts = len(expected_permut_indexes)
    
    print nb_permuts
    
    ################# selecting columns
    
    if len(columns) != 0:
        data_cols = columns
        
    else:
        
        if stop_col == 0:
            data_cols = permut_df.columns[start_col:]
            
        else:
            data_cols = permut_df.columns[start_col:stop_col]
            
    print data_cols
    
    ################## looping over selected columns
    all_p_higher = np.zeros(shape = (len(data_cols)), dtype = 'float64') -1
    all_p_lower = np.zeros(shape = (len(data_cols)), dtype = 'float64') -1
        
    cols = []

    if session_col == -1 or len(permut_df[session_col].unique()) == 1:
        
        for index_col,col in enumerate(data_cols):
            
            print index_col,col
            
            sum_higher = np.sum((permut_df[col].iloc[1:] > permut_df[col].iloc[0]).values.astype(int))
            all_p_higher[index_col] = (sum_higher+1)/float(permut_df[col].shape[0])
            
            sum_lower = np.sum((permut_df[col].iloc[1:] < permut_df[col].iloc[0]).values.astype(int))
            all_p_lower[index_col] = (sum_lower+1)/float(permut_df[col].shape[0])
            
            #print permut_df[col]
            #print permut_df[col].shape[0]
            
            print all_p_higher[index_col] 
            
            cols.append(str(col))
            
            
        print all_p_higher
        print cols
        
        #df_res = pd.DataFrame(all_p_higher.reshape(1,-1),columns = cols)
        #df_res.index = ["Higher"]
            
        #print df_res
        
        #return df_res
    
    else:
        ### all unique values should have 2 different samples
        count_elements = Counter(permut_df[permut_col].values)
        
        print count_elements
        
        if not all(val == 2 for val in count_elements.values()):
            print "Error, all permut indexes should have 2 and only 2 lines: {}".format(count_elements)
            #0/0
            #return pd.DataFrame()
        
        ################################################## computing diff df
            
        
        for index_col,col in enumerate(data_cols):
            
            
            print index_col,col
            
            #print permut_df
            
            print permut_col,session_col
            
        
            df_col = permut_df.pivot(index = permut_col, columns = session_col, values = col)
            
            print df_col
            
            df_col["Diff"] = pd.to_numeric(df_col.iloc[:,0]) - pd.to_numeric(df_col.iloc[:,1])
            
            print df_col["Diff"]
            
            non_nan_indexes, = np.where(np.isnan(df_col["Diff"]) == False)
            
            print non_nan_indexes
            
            diff_col = df_col["Diff"].values[non_nan_indexes]
            
            
            if not all(val == 2 for val in count_elements.values()):
                print "Error, all permut indexes should have 2 and only 2 lines: {}".format(count_elements)
                
                print diff_col
                print diff_col.shape
                
            #assert df_col.shape[0] == 201, "Error with shape {}".format(df_col.shape)
            
            if diff_col.shape[0] == 0:
                
                all_p_higher[index_col] = np.nan
                all_p_lower[index_col] = np.nan
                cols.append(col)
                
                continue
            
            if diff_col[0] > 0:
                sum_higher = np.sum(np.array(diff_col[1:] > diff_col[0],dtype = int))
                print "sum_higher:",sum_higher
                all_p_higher[index_col] = (sum_higher+1)/float(df_col.shape[0])
                
            elif diff_col[0] < 0 :
                sum_lower = np.sum(np.array(diff_col[1:] < diff_col[0],dtype = int))
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
    
def compute_signif_node_prop(orig_df, list_permut_df, columns):

    permut_df = pd.concat(list_permut_df,axis = 0)
    
    print permut_df['Seed']
    
    all_frac_higher = []
    
    for col in columns:
        
        print col
        
        assert col in orig_df.columns, "Error, {} not in orig columns {}".format(col, orig_df.columns)
        assert col in permut_df.columns, "Error, {} not in permut columns {}".format(col, permut_df.columns)
        
                
        def sum_higher(a,b):
            
            print a
            print b
            
            print a[:,None]
            
            def func(el):
                print el[0]
                print b.values
                
                print el[0]<b.values
                
                return np.sum(el[0]<b.values)
            return  np.apply_along_axis(func, 1, a[:,None])

        print orig_df[col]
        print permut_df[col]
        
        frac_higher = np.array(sum_higher(orig_df[col],permut_df[col])+1,dtype = float)/float(len(permut_df.index) +1)

        print frac_higher
        
        all_frac_higher.append(frac_higher)
        
    df_signif = pd.DataFrame(np.transpose(np.array(all_frac_higher)),columns = columns)
    
    return df_signif

##################################### main ###############################################
if __name__ =='__main__':
	
	test1 = isInAlphabeticalOrder(["a","b","c"])
	print test1
	
	test2 = isInAlphabeticalOrder(["ab","ba","ca"])
	print test2
	
