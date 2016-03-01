# -*- coding: utf-8 -*-

import numpy as np
import os

from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec, isdefined
    

from nipype.utils.filemanip import split_filename as split_f
        
############################################################################################### StatsPairBinomial #####################################################################################################

import dmgraphanalysis_nodes.utils_stats as stats
        
class StatsPairBinomialInputSpec(BaseInterfaceInputSpec):
    
    group_coclass_matrix_file1 = File(exists=True,  desc='file of group 1 coclass matrices in npy format', mandatory=True)
    group_coclass_matrix_file2 = File(exists=True,  desc='file of group 2 coclass matrices in npy format', mandatory=True)
    
    conf_interval_binom_fdr = traits.Float(0.05, usedefault = True, desc='Alpha value used as FDR implementation', mandatory=False)
    
class StatsPairBinomialOutputSpec(TraitedSpec):
    
    signif_signed_adj_fdr_mat_file = File(exists=True, desc="int matrix with corresponding codes to significance")
    
class StatsPairBinomial(BaseInterface):
    
    """
    Plot coclassification matrix with igraph
    - labels are optional, 
    - threshold is optional (default, 50 = half the group)
    - coordinates are optional, if no coordiantes are specified, representation in topological (Fruchterman-Reingold) space
    """
    input_spec = StatsPairBinomialInputSpec
    output_spec = StatsPairBinomialOutputSpec

    def _run_interface(self, runtime):
                
        print 'in plot_coclass'
        
        group_coclass_matrix_file1 = self.inputs.group_coclass_matrix_file1
        group_coclass_matrix_file2 = self.inputs.group_coclass_matrix_file2
        conf_interval_binom_fdr = self.inputs.conf_interval_binom_fdr
            

        print "loading group_coclass_matrix1"
        
        group_coclass_matrix1 = np.array(np.load(group_coclass_matrix_file1),dtype = float)
        print group_coclass_matrix1.shape
        
        
        print "loading group_coclass_matrix2"
        
        group_coclass_matrix2 = np.array(np.load(group_coclass_matrix_file2),dtype = float)
        print group_coclass_matrix2.shape
        
        
        print "compute NBS stats"
        
        
        # check input matrices
        Ix,Jx,nx = group_coclass_matrix1.shape
        Iy,Jy,ny = group_coclass_matrix2.shape
        
        assert Ix == Iy
        assert Jx == Jy
        assert Ix == Jx
        assert Iy == Jy
        
        signif_signed_adj_mat  = stats.compute_pairwise_binom_fdr(group_coclass_matrix1,group_coclass_matrix2,conf_interval_binom_fdr)
        
        print 'save pairwise signed stat file'
        
        signif_signed_adj_fdr_mat_file  = os.path.abspath('signif_signed_adj_fdr_'+ str(conf_interval_binom_fdr) +'.npy')
        np.save(signif_signed_adj_fdr_mat_file,signif_signed_adj_mat)
        
        #return signif_signed_adj_fdr_mat_file

            
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["signif_signed_adj_fdr_mat_file"] = os.path.abspath('signif_signed_adj_fdr_'+ str(self.inputs.conf_interval_binom_fdr) +'.npy')
        
        return outputs

############################################################################################### StatsPairTTest #####################################################################################################

import dmgraphanalysis_nodes.utils_stats as stats
        
class StatsPairTTestInputSpec(BaseInterfaceInputSpec):
    
    group_cormat_file1 = File(exists=True,  desc='file of group 1 cormat matrices in npy format', mandatory=True)
    group_cormat_file2 = File(exists=True,  desc='file of group 2 cormat matrices in npy format', mandatory=True)
    
    t_test_thresh_fdr = traits.Float(0.05, usedefault = True, desc='Alpha value used as FDR implementation', mandatory=False)
    
    paired = traits.Bool(True,usedefault = True, desc='Ttest is paired or not', mandatory=False)
    
class StatsPairTTestOutputSpec(TraitedSpec):
    
    signif_signed_adj_fdr_mat_file = File(exists=True, desc="int matrix with corresponding codes to significance")
    
class StatsPairTTest(BaseInterface):
    
    """
    Compute ttest stats between 2 group of matrix 
    - matrix are arranged in group_cormat, with order (Nx,Ny,Nsubj). Nx = Ny (each matricx is square)
    - t_test_thresh_fdr is optional (default, 0.05)
    - paired in indicate if ttest is pairde or not. If paired, both group have the same number of samples
    """
    input_spec = StatsPairTTestInputSpec
    output_spec = StatsPairTTestOutputSpec

    def _run_interface(self, runtime):
                
        print 'in plot_cormat'
        
        group_cormat_file1 = self.inputs.group_cormat_file1
        group_cormat_file2 = self.inputs.group_cormat_file2
        t_test_thresh_fdr = self.inputs.t_test_thresh_fdr
            
        paired = self.inputs.paired
        print "loading group_cormat1"
        
        group_cormat1 = np.array(np.load(group_cormat_file1),dtype = float)
        print group_cormat1.shape
        
        
        print "loading group_cormat2"
        
        group_cormat2 = np.array(np.load(group_cormat_file2),dtype = float)
        print group_cormat2.shape
        
        
        print "compute NBS stats"
        
        
        # check input matrices
        Ix,Jx,nx = group_cormat1.shape
        Iy,Jy,ny = group_cormat2.shape
        
        assert Ix == Iy
        assert Jx == Jy
        assert Ix == Jx
        assert Iy == Jy
        
        signif_signed_adj_mat  = stats.compute_pairwise_ttest_fdr(group_cormat1,group_cormat2,t_test_thresh_fdr,paired)
        
        print 'save pairwise signed stat file'
        
        signif_signed_adj_fdr_mat_file  = os.path.abspath('signif_signed_adj_fdr_'+ str(t_test_thresh_fdr) +'.npy')
        np.save(signif_signed_adj_fdr_mat_file,signif_signed_adj_mat)
        
        #return signif_signed_adj_fdr_mat_file

            
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["signif_signed_adj_fdr_mat_file"] = os.path.abspath('signif_signed_adj_fdr_'+ str(self.inputs.t_test_thresh_fdr) +'.npy')
        
        return outputs

        
############################################################################################### PlotIGraphSignedIntMat #####################################################################################################

from dmgraphanalysis_nodes.plot_igraph import plot_3D_igraph_signed_int_mat
from dmgraphanalysis_nodes.utils import check_np_shapes

#from dmgraphanalysis.utils_plot import plot_cormat
    
    
class PlotIGraphSignedIntMatInputSpec(BaseInterfaceInputSpec):
    
    signed_int_mat_file = File(exists=True,  desc='signed int matrix in npy format', mandatory=True)
    
    labels_file = File(exists=True,  desc='labels of nodes (txt file)', mandatory=False)
    coords_file = File(exists=True,  desc='node coordinates in MNI space (txt file)', mandatory=False)
    
class PlotIGraphSignedIntMatOutputSpec(TraitedSpec):
    
    plot_3D_signed_int_mat_file = File(exists=True, desc="eps file with igraph spatial representation")
    #plot_heatmap_signed_bin_mat_file = File(exists=True, desc="eps file heatmap representation")
    
class PlotIGraphSignedIntMat(BaseInterface):
    
    """
    Plot coclassification matrix with igraph
    - labels are optional, 
    - threshold is optional (default, 50 = half the group)
    - coordinates are optional, if no coordiantes are specified, representation in topological (Fruchterman-Reingold) space
    """
    input_spec = PlotIGraphSignedIntMatInputSpec
    output_spec = PlotIGraphSignedIntMatOutputSpec

    def _run_interface(self, runtime):
                
        print 'in plot_coclass'
        
        signed_int_mat_file = self.inputs.signed_int_mat_file
        labels_file = self.inputs.labels_file
        coords_file = self.inputs.coords_file
            
            
        print 'load bin matrix'
        
        signed_int_mat = np.load(signed_int_mat_file)
        
        print signed_int_mat.shape
        
        if isdefined(labels_file):
            
            print 'loading labels'
            labels = [line.strip() for line in open(labels_file)]
            
        else :
            labels = []
            
        if isdefined(coords_file):
            
            print 'loading coords'
            coords = np.array(np.loadtxt(coords_file),dtype = 'int64')
            
        else :
            coords = np.array([])
            
            
        print coords.shape
        
        
        print 'plotting igraph 3D'
        
        ######## igraph 3D
        plot_3D_signed_int_mat_file = os.path.abspath('plot_igraph_3D_signed_int_mat.eps')
            
        plot_3D_igraph_signed_int_mat(plot_3D_signed_int_mat_file,signed_int_mat,coords,labels)
        
                
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["plot_3D_signed_int_mat_file"] = os.path.abspath('plot_igraph_3D_signed_int_mat.eps')
        #outputs["plot_heatmap_signed_int_mat_file"] = os.path.abspath('heatmap_signed_bin_mat.eps')
        
        return outputs
        
        
        
        
############################################################################################### PrepareCormat #####################################################################################################

from dmgraphanalysis_nodes.utils_cor import return_corres_correl_mat
#,return_hierachical_order
        
class PrepareCormatInputSpec(BaseInterfaceInputSpec):
    
    cor_mat_files = traits.List(File(exists=True), desc='list of all correlation matrice files (in npy format) for each subject', mandatory=True)
    
    coords_files = traits.List(File(exists=True), desc='list of all coordinates in numpy space files (in txt format) for each subject (after removal of non void data)', mandatory=True)
    
    gm_mask_coords_file = File(exists=True, desc='Coordinates in numpy space, corresponding to all possible nodes in the original space', mandatory=True)
    
    
class PrepareCormatOutputSpec(TraitedSpec):
    
    group_cormat_file = File(exists=True, desc="all cormat matrices of the group in .npy (pickle format)")
    
    avg_cormat_file = File(exists=True, desc="average of cormat matrix of the group in .npy (pickle format)")
    
    group_vect_file = File(exists=True, desc="degree (?) by nodes * indiv of the group in .npy (pickle format)")
    
    
class PrepareCormat(BaseInterface):
    
    """
    Extract mean time series from a labelled mask in Nifti Format where the voxels of interest have values 1
    """
    input_spec = PrepareCormatInputSpec
    output_spec = PrepareCormatOutputSpec

    def _run_interface(self, runtime):
                
        print 'in prepare_coclass'
        cor_mat_files = self.inputs.cor_mat_files
        coords_files = self.inputs.coords_files
        gm_mask_coords_file = self.inputs.gm_mask_coords_file
        
        
    #import numpy as np
    #import os

    ##import nibabel as nib
        print 'loading gm mask corres'
        
        gm_mask_coords = np.loadtxt(gm_mask_coords_file)
        
        print gm_mask_coords.shape
            
        #### read matrix from the first group
        #print Z_cor_mat_files
        
        sum_cormat = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = float)
        print sum_cormat.shape
        
                
        group_cormat = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0],len(cor_mat_files)),dtype = float)
        print group_cormat.shape
        
        
        group_vect = np.zeros((gm_mask_coords.shape[0],len(cor_mat_files)),dtype = float)
        print group_vect.shape
        
        if len(cor_mat_files) != len(coords_files):
            print "warning, length of cor_mat_files, coords_files are imcompatible {} {} {}".format(len(cor_mat_files),len(coords_files))
        
        for index_file in range(len(cor_mat_files)):
            
            print cor_mat_files[index_file]
            
            if os.path.exists(cor_mat_files[index_file]) and os.path.exists(coords_files[index_file]):
            
                Z_cor_mat = np.load(cor_mat_files[index_file])
                print Z_cor_mat.shape
                
                
                coords = np.loadtxt(coords_files[index_file])
                print coords.shape
                
                
                
                corres_cor_mat,possible_edge_mat = return_corres_correl_mat(Z_cor_mat,coords,gm_mask_coords)
                
                print corres_cor_mat.shape
                print group_cormat.shape
                
                sum_cormat += corres_cor_mat
                
                group_cormat[:,:,index_file] = corres_cor_mat
                
                group_vect[:,index_file] = np.sum(corres_cor_mat,axis = 0)
                
                
            else:
                print "Warning, one or more files between " + cor_mat_files[index_file] + ', ' + coords_files[index_file] + " do not exists"
            
            
        group_cormat_file= os.path.abspath('group_cormat.npy')
        
        np.save(group_cormat_file,group_cormat)
        
            
        group_vect_file= os.path.abspath('group_vect.npy')
        
        np.save(group_vect_file,group_vect)
        
            
        print 'saving cor_mat matrix'
        
        avg_cormat_file = os.path.abspath('avg_cormat.npy')
        
        if (len(cor_mat_files) != 0):
        
                avg_cormat = sum_cormat /len(cor_mat_files)
                
                np.save(avg_cormat_file,avg_cormat)
        
        
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["group_cormat_file"] =os.path.abspath('group_cormat.npy')
        
        outputs["avg_cormat_file"] = os.path.abspath('avg_cormat.npy')
        
        outputs["group_vect_file"] = os.path.abspath('group_vect.npy')
        
        return outputs