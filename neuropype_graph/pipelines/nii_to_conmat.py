# -*- coding: utf-8 -*-
"""
From nifti file to conmat
"""
import sys
import time

import numpy as np
import scipy.sparse as sp 

import nipype.pipeline.engine as pe

import nipype.interfaces.utility  as niu
import nipype.interfaces.spm.utils as spmu

from neuropype_graph.nodes.correl_mat import IntersectMask,ExtractTS,ExtractMeanTS,RegressCovar,FindSPMRegressor,MergeRuns,ComputeConfCorMat

from neuropype_graph.utils import show_files,show_length

#import imp

#try:
    #imp.find_module('igraph')
    #can_plot_igraph = True
    #from neuropype_graph.nodes.igraph_plots import PlotIGraphModules

#except ImportError:
    #can_plot_igraph = False
    

def create_pipeline_nii_to_conmat(main_path,ROI_mask_file,ROI_coords_file,ROI_MNI_coords_file,ROI_labels_file, pipeline_name = "nii_to_conmat",conf_interval_prob = 0.05):


    pipeline = pe.Workflow(name=pipeline_name)
    pipeline.base_dir = main_path
    
    inputnode = pe.Node(niu.IdentityInterface(fields=['nii_4D_file','rp_file','gm_anat_file','wm_anat_file','csf_anat_file']),
                        name='inputnode')
     
    ###### Preprocess pipeline,
    filter_ROI_mask_with_GM = pe.Node(interface = IntersectMask(),name = 'filter_ROI_mask_with_GM')
    
    filter_ROI_mask_with_GM.inputs.indexed_rois_file = ROI_mask_file
    
    filter_ROI_mask_with_GM.inputs.coords_rois_file = ROI_coords_file
    filter_ROI_mask_with_GM.inputs.MNI_coords_rois_file = ROI_MNI_coords_file    
    filter_ROI_mask_with_GM.inputs.labels_rois_file = ROI_labels_file
    
    
    pipeline.connect(inputnode, 'gm_anat_file', filter_ROI_mask_with_GM, 'filter_mask_file')
    
    #### Nodes version: use min_BOLD_intensity and return coords where signal is strong enough 
    extract_mean_ROI_ts = pe.Node(interface = ExtractTS(plot_fig = False),name = 'extract_mean_ROI_ts')
    
    
    #extract_mean_ROI_ts.inputs.indexed_rois_file = ROI_mask_file
    #extract_mean_ROI_ts.inputs.coord_rois_file = ROI_coords_file
    #extract_mean_ROI_ts.inputs.min_BOLD_intensity = min_BOLD_intensity
    
    pipeline.connect(inputnode,'nii_4D_file', extract_mean_ROI_ts, 'file_4D')
    pipeline.connect(filter_ROI_mask_with_GM, 'filtered_indexed_rois_file', extract_mean_ROI_ts, 'indexed_rois_file')
    pipeline.connect(filter_ROI_mask_with_GM, 'filtered_coords_rois_file', extract_mean_ROI_ts, 'coord_rois_file')
    
    
    #### reslice white_matter_signal
    reslice_wm = pe.Node(interface = spmu.Reslice(), name = 'reslice_wm')    
    reslice_wm.inputs.space_defining = ROI_mask_file
    
    pipeline.connect(inputnode, 'wm_anat_file', reslice_wm, 'in_file')
    
    #### extract white matter signal
    compute_wm_ts = pe.Node(interface = ExtractMeanTS(plot_fig = False),name = 'extract_wm_ts')
    compute_wm_ts.inputs.suffix = 'wm'
    
    pipeline.connect(inputnode,'nii_4D_file', compute_wm_ts, 'file_4D')
    pipeline.connect(reslice_wm, 'out_file', compute_wm_ts, 'filter_mask_file')
    
    
    #### reslice csf
    reslice_csf = pe.Node(interface = spmu.Reslice(), name = 'reslice_csf')    
    reslice_csf.inputs.space_defining = ROI_mask_file
    
    pipeline.connect(inputnode, 'csf_anat_file', reslice_csf, 'in_file')
    
    
    #### extract csf signal
    compute_csf_ts = pe.Node(interface = ExtractMeanTS(plot_fig = False),name = 'extract_csf_ts')
    compute_csf_ts.inputs.suffix = 'csf'
    
    pipeline.connect(inputnode,'nii_4D_file', compute_csf_ts, 'file_4D')
    pipeline.connect(reslice_csf, 'out_file', compute_csf_ts, 'filter_mask_file')
    
    
    
    #### regress covariates
    
    ### use R linear model to regress movement parameters, white matter and ventricule signals, and compute Z-score of the residuals
    #regress_covar = pe.MapNode(interface = RegressCovar(filtered = False, normalized = False),iterfield = ['masked_ts_file','rp_file','mean_wm_ts_file','mean_csf_ts_file'],name='regress_covar')
    regress_covar = pe.Node(interface = RegressCovar(),iterfield = ['masked_ts_file','rp_file','mean_wm_ts_file','mean_csf_ts_file'],name='regress_covar')
    
    pipeline.connect(extract_mean_ROI_ts, ('mean_masked_ts_file',show_files), regress_covar, 'masked_ts_file')
    pipeline.connect(inputnode, 'rp_file', regress_covar, 'rp_file')

    pipeline.connect(compute_wm_ts, 'mean_masked_ts_file', regress_covar, 'mean_wm_ts_file')
    pipeline.connect(compute_csf_ts, 'mean_masked_ts_file', regress_covar, 'mean_csf_ts_file')
    
    
    
    
    #### extract regressor of interest from SPM.mat
    #extract_cond = pe.Node(interface = FindSPMRegressor(only_positive_values = True),name='extract_cond')
    
    #pipeline.connect(inputnode, ('spm_mat_file',show_files), extract_cond, 'spm_mat_file')
    #pipeline.connect(infosource, 'cond', extract_cond, 'regressor_name')
    
    ##extract_cond.inputs.run_index = funct_run_indexs
    
    
    ##### merge_runs new version: merge also coords (if different between sessions)
    ##merge_runs = pe.Node(interface = MergeRuns(),name='merge_runs')
    
    ##pipeline.connect(extract_mean_ROI_ts, ('subj_coord_rois_file',show_length), merge_runs, 'coord_rois_files')
    ##pipeline.connect(regress_covar, ('resid_ts_file',show_length), merge_runs, 'ts_files')
    ##pipeline.connect(extract_cond, ('regressor_file',show_length), merge_runs, 'regressor_files')
    
    ##################################### compute correlations ####################################################
    
    ########### confidence interval 
    
    compute_conf_cor_mat = pe.Node(interface = ComputeConfCorMat(),name='compute_conf_cor_mat')
    
    compute_conf_cor_mat.inputs.conf_interval_prob = conf_interval_prob
    
    pipeline.connect(regress_covar, ('resid_ts_file',show_length), compute_conf_cor_mat, 'ts_file')
    
    
    ###plot_hist_conf = pe.Node(Function(input_names=['cor_mat_file','Z_cor_mat_file','conf_cor_mat_file'],output_names=['plot_hist_cor_mat_file','plot_heatmap_cor_mat_file','plot_hist_cor_mat_file','plot_heatmap_cor_mat_file','plot_hist_conf_cor_mat_file','plot_heatmap_conf_cor_mat_file'],function=plot_hist_conf_cor_mat),name='plot_hist_conf')
    
    ###pipeline.connect(compute_conf_cor_mat, 'cor_mat_file',plot_hist_conf,'cor_mat_file')
    ###pipeline.connect(compute_conf_cor_mat, 'Z_cor_mat_file',plot_hist_conf,'Z_cor_mat_file')
    ###pipeline.connect(compute_conf_cor_mat, 'conf_cor_mat_file',plot_hist_conf,'conf_cor_mat_file')
    
    return pipeline
