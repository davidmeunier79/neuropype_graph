# -*- coding: utf-8 -*-
"""
Support function for net handling
"""
import sys
import time

import numpy as np
import scipy.sparse as sp 

import nipype.pipeline.engine as pe
    
from dmgraphanalysis_nodes.nodes.modularity import ComputeIntNetList

from dmgraphanalysis_nodes.nodes.modularity import PrepRada,CommRada,ComputeNodeRoles
from dmgraphanalysis_nodes.nodes.modularity import NetPropRada
 
import imp

try:
    imp.find_module('igraph')
    can_plot_igraph = True
    from dmgraphanalysis_nodes.nodes.igraph_plots import PlotIGraphModules

except ImportError:
    can_plot_igraph = False
    
    
 
#def create_pipeline_conmat_to_graph_density(correl_analysis_name,main_path,radatools_path,con_den = 1.0,multi = False,mod = True):

    #pipeline = pe.Workflow(name=correl_analysis_name)
    #pipeline.base_dir = main_path
    
    #if multi == False:
        
        ################################################# density-based graphs
        
        ##### net_list
        #compute_net_List = pe.Node(interface = ComputeNetList(),name='compute_net_List')
        #compute_net_List.inputs.density = con_den
        
        ##pipeline.connect(convert_mat, 'conmat_file',compute_net_List, 'Z_cor_mat_file')
        
        ###### radatools ################################################################

        #### prepare net_list for radatools processing  
        #prep_rada = pe.Node(interface = PrepRada(),name='prep_rada',iterfield = ["net_List_file"])
        #prep_rada.inputs.radatools_path = radatools_path
        
        #pipeline.connect(compute_net_List, 'net_List_file', prep_rada, 'net_List_file')
    
        #if mod == True:
                
            #### compute community with radatools
            #community_rada = pe.Node(interface = CommRada(), name='community_rada')
            ##community_rada.inputs.optim_seq = radatools_optim
            #community_rada.inputs.radatools_path = radatools_path
            
            #pipeline.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
            
            
            
            #### node roles
            #node_roles = pe.Node(interface = ComputeNodeRoles(role_type = "4roles"), name='node_roles')
            
            #pipeline.connect( prep_rada, 'Pajek_net_file',node_roles,'Pajek_net_file')
            #pipeline.connect( community_rada, 'rada_lol_file',node_roles,'rada_lol_file')
            
            ##### plot_igraph_modules_rada
            #plot_igraph_modules_rada = pe.Node(interface = PlotIGraphModules(),name='plot_igraph_modules_rada')
            
            #pipeline.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada,'Pajek_net_file')
            #pipeline.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada,'rada_lol_file')
            
            #pipeline.connect(node_roles, 'node_roles_file',plot_igraph_modules_rada,'node_roles_file')
    
        ############# compute network properties with rada
        #net_prop = pe.Node(interface = NetPropRada(optim_seq = "A"), name = 'net_prop')
        #net_prop.inputs.radatools_path = radatools_path
        
        #pipeline.connect(prep_rada, 'Pajek_net_file',net_prop,'Pajek_net_file')
        
        
    #else:
        
        ################################################# density-based graphs #################################################
        
        ##### net_list
        #compute_net_List = pe.MapNode(interface = ComputeNetList(),name='compute_net_List',iterfield = ["Z_cor_mat_file"])
        #compute_net_List.inputs.density = con_den
        
        ##pipeline.connect(convert_mat, 'conmat_file',compute_net_List, 'Z_cor_mat_file')
        
        ###### radatools ################################################################

        #### prepare net_list for radatools processing  
        #prep_rada = pe.MapNode(interface = PrepRada(),name='prep_rada',iterfield = ["net_List_file"])
        #prep_rada.inputs.radatools_path = radatools_path
        
        #pipeline.connect(compute_net_List, 'net_List_file', prep_rada, 'net_List_file')
        
        
        #if mod == True:
                
            #### compute community with radatools
            #community_rada = pe.MapNode(interface = CommRada(), name='community_rada',iterfield = ["Pajek_net_file"])
            ##community_rada.inputs.optim_seq = radatools_optim
            #community_rada.inputs.radatools_path = radatools_path
            
            #pipeline.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
            
            #### node roles
            #node_roles = pe.Node(interface = ComputeNodeRoles(role_type = "4roles"), name='node_roles')
            
            #pipeline.connect( prep_rada, 'Pajek_net_file',node_roles,'Pajek_net_file')
            #pipeline.connect( community_rada, 'rada_lol_file',node_roles,'rada_lol_file')
            
            ##### plot_igraph_modules_rada
            #plot_igraph_modules_rada = pe.MapNode(interface = PlotIGraphModules(),name='plot_igraph_modules_rada',iterfield = ['Pajek_net_file','rada_lol_file','node_roles_file'])
            
            #pipeline.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada,'Pajek_net_file')
            #pipeline.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada,'rada_lol_file')
            
            #pipeline.connect(node_roles, 'node_roles_file',plot_igraph_modules_rada,'node_roles_file')
            
    
        ############# compute network properties with rada
        #net_prop = pe.MapNode(interface = NetPropRada(optim_seq = "A"), name = 'net_prop',iterfield = ["Pajek_net_file"])
        #net_prop.inputs.radatools_path = radatools_path
        
        #pipeline.connect(prep_rada, 'Pajek_net_file',net_prop,'Pajek_net_file')


    #return pipeline

################################################ threshold-based graphs
    
def create_pipeline_intmat_to_graph_threshold(analysis_name,main_path,radatools_path,threshold = 50, plot = False):

    pipeline = pe.Workflow(name=analysis_name)
    pipeline.base_dir = main_path


    if plot==True and can_plot_igraph==False:
        
        plot = False
        
        
    ### compute Z_list from coclass matrix
    compute_list_norm_coclass = pe.Node(interface = ComputeIntNetList(),name='compute_list_norm_coclass')
    compute_list_norm_coclass.inputs.threshold = threshold
    
    ############################################### radatools ################################################################
    
    ###--- prepare net_list for radatools processing
    prep_rada = pe.Node(interface = PrepRada(),name='prep_rada')
    prep_rada.inputs.radatools_path = radatools_path
    
    pipeline.connect(compute_list_norm_coclass, 'net_List_file', prep_rada, 'net_List_file')
    
    
    ### compute community with radatools
    community_rada = pe.Node(interface = CommRada(), name='community_rada',iterfield = ["Pajek_net_file"])
    #community_rada.inputs.optim_seq = radatools_optim
    community_rada.inputs.radatools_path = radatools_path
    
    pipeline.connect( prep_rada, 'Pajek_net_file',community_rada,'Pajek_net_file')
    

    
    ### node roles
    node_roles = pe.Node(interface = ComputeNodeRoles(role_type = "4roles"), name='node_roles')
    
    pipeline.connect( prep_rada, 'Pajek_net_file',node_roles,'Pajek_net_file')
    pipeline.connect( community_rada, 'rada_lol_file',node_roles,'rada_lol_file')
    
    if plot == True:
            
        #### plot_igraph_modules_rada
        plot_igraph_modules_rada = pe.Node(interface = PlotIGraphModules(),name='plot_igraph_modules_rada')
        
        pipeline.connect(prep_rada, 'Pajek_net_file',plot_igraph_modules_rada,'Pajek_net_file')
        pipeline.connect(community_rada, 'rada_lol_file',plot_igraph_modules_rada,'rada_lol_file')
        
        pipeline.connect(node_roles, 'node_roles_file',plot_igraph_modules_rada,'node_roles_file')


    ############################################ compute network properties with rada ############################################
    net_prop = pe.Node(interface = NetPropRada(optim_seq = "A"), name = 'net_prop')
    net_prop.inputs.radatools_path = radatools_path
    
    pipeline.connect(prep_rada, 'Pajek_net_file',net_prop,'Pajek_net_file')
    
    return pipeline

   