# -*- coding: utf-8 -*-

"""
Definition of Nodes for computing correlation matrices 
"""

#import nipy.labs.statistical_mapping as stat_map

#import itertools as iter
    
#import scipy.spatial.distance as dist
    
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec, isdefined
    
from nipype.utils.filemanip import split_filename as split_f
    
import nibabel as nb
import numpy as np
import os

import nibabel as nib

from neuropype_graph.utils_plot import plot_signals,plot_sep_signals
        
######################################################################################## ExtractTS ##################################################################################################################

from neuropype_graph.utils_cor import mean_select_indexed_mask_data

class ExtractTSInputSpec(BaseInterfaceInputSpec):
    indexed_rois_file = File(exists=True, desc='indexed mask where all voxels belonging to the same ROI have the same value (! starting from 1)', mandatory=True)
    
    file_4D = File(exists=True, desc='4D volume to be extracted', mandatory=True)
    
    coord_rois_file = File(desc='ROI coordinates')
    
    min_BOLD_intensity = traits.Float(50.0, desc='BOLD signal below this value will be set to zero',usedefault = True)

    percent_signal = traits.Float(0.5, desc  = "Percent of voxels in a ROI with signal higher that min_BOLD_intensity to keep this ROI",usedefault = True)
    
    plot_fig = traits.Bool(False, desc = "Plotting mean signal or not", usedefault = True)
    
    
class ExtractTSOutputSpec(TraitedSpec):
    
    mean_masked_ts_file = File(exists=True, desc="mean ts in .npy (pickle format)")
    
    subj_coord_rois_file = File(exists=True, desc="ROI coordinates retained for the subject")
    

class ExtractTS(BaseInterface):
    
    """Extract time series from a labelled mask in Nifti Format where all ROIs have the same index"""

    input_spec = ExtractTSInputSpec
    output_spec = ExtractTSOutputSpec

    def _run_interface(self, runtime):
            
        #import os
        #import numpy as np
        #import nibabel as nib
                
        coord_rois_file = self.inputs.coord_rois_file
        indexed_rois_file = self.inputs.indexed_rois_file
        file_4D = self.inputs.file_4D
        min_BOLD_intensity = self.inputs.min_BOLD_intensity
        
        plot_fig = self.inputs.plot_fig
        
        ## loading ROI coordinates
        coord_rois = np.loadtxt(coord_rois_file)
        
        print "coord_rois: " 
        print coord_rois.shape
        
        ## loading ROI indexed mask
        indexed_rois_img = nib.load(indexed_rois_file)
        
        indexed_mask_rois_data = indexed_rois_img.get_data()
        
        #print "indexed_mask_rois_data: "
        #print indexed_mask_rois_data.shape
        
        ### loading time series
        orig_ts = nib.load(file_4D).get_data()
        
        print "orig_ts shape:"
        print orig_ts.shape
            
        mean_masked_ts,subj_coord_rois = mean_select_indexed_mask_data(orig_ts,indexed_mask_rois_data,coord_rois,min_BOLD_intensity, percent_signal = 0.5)
        
        mean_masked_ts = np.array(mean_masked_ts,dtype = 'f')
        subj_coord_rois = np.array(subj_coord_rois,dtype = 'float')
        
        print mean_masked_ts.shape
            
        ### saving time series
        mean_masked_ts_file = os.path.abspath("mean_masked_ts.txt")
        np.savetxt(mean_masked_ts_file,mean_masked_ts,fmt = '%.3f')
        
        ### saving subject ROIs
        subj_coord_rois_file = os.path.abspath("subj_coord_rois.txt")
        np.savetxt(subj_coord_rois_file,subj_coord_rois,fmt = '%.3f')
        
        if plot_fig == True:
                
            print "plotting mean_masked_ts"
            
            plot_mean_masked_ts_file = os.path.abspath('mean_masked_ts.eps')    
            
            plot_signals(plot_mean_masked_ts_file,mean_masked_ts)
            
        return runtime
        
        #return mean_masked_ts_file,subj_coord_rois_file
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["mean_masked_ts_file"] = os.path.abspath("mean_masked_ts.txt")
        outputs["subj_coord_rois_file"] = os.path.abspath("subj_coord_rois.txt")
    
        return outputs


######################################################################################## IntersectMask ##################################################################################################################

from neuropype_graph.utils_cor import mean_select_indexed_mask_data

class IntersectMaskInputSpec(BaseInterfaceInputSpec):
    
    
    indexed_rois_file = File(exists=True, desc='nii file with indexed mask where all voxels belonging to the same ROI have the same value (! starting from 1)', mandatory=True)
    filter_mask_file = File( exists=True, desc='nii file with (binary) mask - e.g. grey matter mask', mandatory=True)
    
    coords_rois_file = File(desc='ijk coords txt file')
    labels_rois_file = File( desc='labels txt file')
    MNI_coords_rois_file = File(desc='MNI coords txt file')
    
    filter_thr = traits.Float(0.99, usedefault = True, desc='Value to threshold filter_mask')
                          
class IntersectMaskOutputSpec(TraitedSpec):
    
    filtered_indexed_rois_file = File(exists=True, desc='nii file with indexed mask where all voxels belonging to the same ROI have the same value (! starting from 1)')
    
    
    filtered_coords_rois_file = File(exists=False, desc='filtered ijk coords txt file')    
    #
    #filtered_labels_rois_file = File(exists=False, desc='filtered labels txt file')
    #filtered_MNI_coords_rois_file = File(exists=False, desc='filtered MNI coords txt file')
    
    


class IntersectMask(BaseInterface):
    
    """
    Keep only values of indexed mask where filtermask is present
    Optionnally, keep only ijk_coords, MNI_coords and labels that are kept in filtered mask
    """

    input_spec = IntersectMaskInputSpec
    output_spec = IntersectMaskOutputSpec

    def _run_interface(self, runtime):
            
        #import os
        #import numpy as np
        import nibabel as nib
        
        import nipype.interfaces.spm as spm

        #from neuropype_graph.utils_plot import plot_signals
        
        indexed_rois_file = self.inputs.indexed_rois_file
        filter_mask_file = self.inputs.filter_mask_file
        coords_rois_file = self.inputs.coords_rois_file
        labels_rois_file = self.inputs.labels_rois_file
        MNI_coords_rois_file = self.inputs.MNI_coords_rois_file
        filter_thr = self.inputs.filter_thr 
        
        print filter_thr
        
        ## loading ROI indexed mask
        indexed_rois_img = nib.load(indexed_rois_file)
        
        indexed_rois_data = indexed_rois_img.get_data()
        
        print "indexed_rois_data: "
        print indexed_rois_data.shape
        
        ### loading time series
        filter_mask_data = nib.load(filter_mask_file).get_data()
        
        print "filter_mask_data shape:"
        print filter_mask_data.shape
            
        if filter_mask_data.shape != indexed_rois_data.shape:
            
            print "reslicing filtered_mask"
            
            
            reslice_filter_mask = spm.Reslice()
            reslice_filter_mask.inputs.in_file = filter_mask_file
            reslice_filter_mask.inputs.space_defining = indexed_rois_file
            #reslice_filter_mask.inputs.out_file = os.path.abspath("resliced_filter_mask.nii")
            
            resliced_filter_mask_file =  reslice_filter_mask.run().outputs.out_file

            filter_mask_data = nib.load(resliced_filter_mask_file).get_data()
        
        print filter_mask_data.shape
        
        print np.unique(filter_mask_data)
        
        print np.sum(filter_mask_data > filter_thr)
        
        filter_mask_data[filter_mask_data > filter_thr] = 1.0
        filter_mask_data[filter_mask_data <= filter_thr] = 0.0
        
        #print np.unique(filter_mask_data)
        print "indexed_rois_data:"
        print np.unique(indexed_rois_data)        
        print len(np.unique(indexed_rois_data))
        
        filtered_indexed_rois_data = np.array(filter_mask_data * (indexed_rois_data.copy()+1) -1,dtype = 'int64')
        
        print "filtered_indexed_rois_data:"
        print np.unique(filtered_indexed_rois_data)
        print len(np.unique(filtered_indexed_rois_data))
            
        filtered_indexed_rois_img_file = os.path.abspath("filtered_indexed_rois.nii")
        nib.save(nib.Nifti1Image(filtered_indexed_rois_data,indexed_rois_img.get_affine(),indexed_rois_img.get_header()),filtered_indexed_rois_img_file)
    
        print "index_corres:"
        #index_corres = np.unique(filtered_indexed_rois_data)[:-1]
        index_corres = np.unique(filtered_indexed_rois_data)[1:]
        
        print index_corres
        print len(index_corres)
        
        print "reorder_indexed_rois:"
        reorder_indexed_rois_data = np.zeros(shape = filtered_indexed_rois_data.shape,dtype = 'int64') - 1
        
        for i,index in enumerate(index_corres):
            
            print i,index
            
            if np.sum(np.array(filtered_indexed_rois_data == index,dtype = int)) != 0:
                
                print np.sum(np.array(filtered_indexed_rois_data == index,dtype = int))
                reorder_indexed_rois_data[filtered_indexed_rois_data == index] = i
            else:
                print "Warning could not find value %d in filtered_indexed_rois_data"%index
            
                                       
        print np.unique(reorder_indexed_rois_data)  
        
        reorder_indexed_rois_img_file = os.path.abspath("reorder_filtered_indexed_rois.nii")
        nib.save(nib.Nifti1Image(reorder_indexed_rois_data,indexed_rois_img.get_affine(),indexed_rois_img.get_header()),reorder_indexed_rois_img_file)
    
    
        if isdefined(coords_rois_file):
            
            ## loading ROI coordinates
            print "coords_rois_file: "
            print coords_rois_file
            coords_rois = np.loadtxt(coords_rois_file)
            
            print "coords_rois: " 
            print coords_rois.shape
            
            filtered_coords_rois = coords_rois[index_corres,:]
            
            print "filtered_coords_rois: " 
            print filtered_coords_rois
            
            filtered_coords_rois_file = os.path.abspath("filtered_coords_rois.txt")
            np.savetxt(filtered_coords_rois_file,filtered_coords_rois, fmt = "%d")
            
        if isdefined(MNI_coords_rois_file):
            
            ## loading ROI coordinates
            MNI_coords_rois = np.loadtxt(MNI_coords_rois_file)
            
            print "MNI_coords_rois: " 
            print MNI_coords_rois.shape
            
            filtered_MNI_coords_rois = MNI_coords_rois[index_corres,:]
            
            print "filtered_MNI_coords_rois: " 
            print filtered_MNI_coords_rois
            
            filtered_MNI_coords_rois_file = os.path.abspath("filtered_MNI_coords_rois.txt")
            np.savetxt(filtered_MNI_coords_rois_file,filtered_MNI_coords_rois, fmt = "%f")
            
        if isdefined(labels_rois_file):
    
    
            print 'extracting node labels'
                
            labels_rois = [line.strip() for line in open(labels_rois_file)]
            print labels_rois
            print len(labels_rois)
            
            np_labels_rois = np.array(labels_rois,dtype = 'str')
            
                
            filtered_labels_rois = np_labels_rois[index_corres]
            
            print "filtered_coords_rois: " 
            print filtered_coords_rois
            
            filtered_labels_rois_file = os.path.abspath("filtered_labels_rois.txt")
            np.savetxt(filtered_labels_rois_file,filtered_labels_rois, fmt = "%s")
            
        return runtime
        
        #return mean_masked_ts_file,subj_coord_rois_file
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["filtered_indexed_rois_file"] = os.path.abspath("reorder_filtered_indexed_rois.nii")
    
        if isdefined(self.inputs.coords_rois_file):
            outputs["filtered_coords_rois_file"] = os.path.abspath("filtered_coords_rois.txt")
            
        #outputs["filtered_MNI_coords_rois_file"] = os.path.abspath("filtered_MNI_coords_rois.txt")
        #outputs["filtered_labels_rois_file"] = os.path.abspath("filtered_labels_rois.txt")
    
        return outputs



############################################################################################### ExtractMeanTS #####################################################################################################

from neuropype_graph.utils_cor import mean_select_mask_data

class ExtractMeanTSInputSpec(BaseInterfaceInputSpec):
    mask_file = File(xor = ['filter_mask_file'], exists=True, desc='mask file where all voxels belonging to the selected region have index 1', mandatory=True)
    
    filter_mask_file = File(xor = ['mask_file'],requires = ['filter_thr'], exists=True, desc='mask file where all voxels belonging to the selected region have values higher than threshold', mandatory=True)
    
    filter_thr = traits.Float(0.99, usedefault = True, desc='Value to threshold filter_mask')
    
    file_4D = File(exists=True, desc='4D volume to be extracted', mandatory=True)
    
    suffix = traits.String(desc='Suffix added to describe the extracted time series',mandatory=False)

    plot_fig = traits.Bool(False, desc = "Plotting mean signal or not", usedefault = True)
    
class ExtractMeanTSOutputSpec(TraitedSpec):
    
    mean_masked_ts_file = File(exists=True, desc="mean ts in .npy (pickle format)")
    

class ExtractMeanTS(BaseInterface):
    
    """
    Extract mean time series from a labelled mask in Nifti Format where the voxels of interest have values 1
    """
    input_spec = ExtractMeanTSInputSpec
    output_spec = ExtractMeanTSOutputSpec

    def _run_interface(self, runtime):
                
        print 'in select_ts_with_mask'
        
        file_4D = self.inputs.file_4D
        mask_file = self.inputs.mask_file
        filter_mask_file = self.inputs.filter_mask_file
        filter_thr = self.inputs.filter_thr
        plot_fig = self.inputs.plot_fig
        
        
        if isdefined(self.inputs.suffix):
            suffix = self.inputs.suffix
        else:
            suffix = "suf"
            
            
        print "loading img data " + file_4D

        ### Reading 4D volume file to extract time series
        img = nib.load(file_4D)
        img_data = img.get_data()
        
        print img_data.shape

        
        ### Reading 4D volume file to extract time series
        if isdefined(mask_file):
        
            print "loading mask data " + mask_file

            mask_data = nib.load(mask_file).get_data()
            
        elif isdefined(filter_mask_file):
            print "loading filter mask data " + filter_mask_file

            filter_mask_data = nib.load(filter_mask_file).get_data()
            mask_data = np.zeros(shape = filter_mask_data.shape, dtype = 'int')
            
            mask_data[filter_mask_data > filter_thr] = 1
            
        print np.unique(mask_data)
        print mask_data.shape
        
        print "mean_select_mask_data"
        
        ### Retaining only time series who are within the mask + non_zero
        mean_masked_ts = mean_select_mask_data(img_data,mask_data)
        
        print mean_masked_ts
        print mean_masked_ts.shape
        
        print "saving mean_masked_ts"
        mean_masked_ts_file = os.path.abspath('mean_' + suffix + '_ts.txt')    
        np.savetxt(mean_masked_ts_file,mean_masked_ts,fmt = '%.3f')
        
        if plot_fig == True:
                
            print "plotting mean_masked_ts"
            
            plot_mean_masked_ts_file = os.path.abspath('mean_' + suffix + '_ts.eps')    
            
            plot_signals(plot_mean_masked_ts_file,mean_masked_ts)
            
        return runtime
        
        #return mean_masked_ts_file,subj_coord_rois_file
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        if isdefined(self.inputs.suffix):
        
            suffix = self.inputs.suffix
        
        else:
            
            suffix = "suf"
            
        outputs["mean_masked_ts_file"] = os.path.abspath('mean_' + suffix + '_ts.txt')
    
        return outputs
        
        
        
######################################################################################## ConcatTS ##################################################################################################################

class ConcatTSInputSpec(BaseInterfaceInputSpec):
    
    all_ts_file = File(exists=True, desc='npy file containing all ts to be concatenated', mandatory=True)
    
class ConcatTSOutputSpec(TraitedSpec):
    
    concatenated_ts_file = File(exists=True, desc="ts after concatenation")
        

class ConcatTS(BaseInterface):
    
    """Concenate time series """

    input_spec = ConcatTSInputSpec
    output_spec = ConcatTSOutputSpec

    def _run_interface(self, runtime):
            
        #import os
        #import numpy as np
        #import nibabel as nib
        
        #from neuropype_graph.utils_plot import plot_signals
        
        all_ts_file = self.inputs.all_ts_file
        
        
        ## loading ROI coordinates
        all_ts = np.load(all_ts_file)
        
        print "all_ts: " 
        print all_ts.shape
        
        concatenated_ts = all_ts.swapaxes(1,0).reshape(all_ts.shape[1],-1)
        
        print concatenated_ts.shape
        
        ### saving time series
        concatenated_ts_file = os.path.abspath("concatenated_ts.npy")
        np.save(concatenated_ts_file,concatenated_ts)
        
        return runtime
        
        #return mean_masked_ts_file,subj_coord_rois_file
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["concatenated_ts_file"] = os.path.abspath("concatenated_ts.npy")
    
        return outputs

######################################################################################## MergeTS ##################################################################################################################

class MergeTSInputSpec(BaseInterfaceInputSpec):
    
    all_ts_files = traits.List(File(exists=True), desc='list of npy files containing all ts to be merged', mandatory=True)
             
class MergeTSOutputSpec(TraitedSpec):
    
    merged_ts_file = File(exists=True, desc="ts after merge")
        

class MergeTS(BaseInterface):
    
    """Merges time series from several files """

    input_spec = MergeTSInputSpec
    output_spec = MergeTSOutputSpec

    def _run_interface(self, runtime):
            
        all_ts_files = self.inputs.all_ts_files
        
        print all_ts_files

        for i,all_ts_file in enumerate(all_ts_files):
        
            all_ts = np.load(all_ts_file)
        
            concatenated_ts = all_ts.swapaxes(1,0).reshape(all_ts.shape[1],-1)
        
            print concatenated_ts.shape

            if len(concatenated_ts.shape) > 1:

                if i == 0:
                    merged_ts = concatenated_ts.copy()
                    print merged_ts.shape
                else:
                    merged_ts = np.concatenate((merged_ts,concatenated_ts),axis = 1)
                    print merged_ts.shape

        merged_ts_file = os.path.abspath("merged_ts.npy")
        np.save(merged_ts_file,merged_ts)
        
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["merged_ts_file"] = os.path.abspath("merged_ts.npy")
    
        return outputs

        
######################################################################################## SeparateTS ##################################################################################################################

class SeparateTSInputSpec(BaseInterfaceInputSpec):
    
    all_ts_file = File(exists=True, desc='npy file containing all ts to be concatenated', mandatory=True)
    
class SeparateTSOutputSpec(TraitedSpec):
    
    separated_ts_files = traits.List(File(exists=True), desc="ts files after separation")
    
class SeparateTS(BaseInterface):
    
    """Extract time series from a labelled mask in Nifti Format where all ROIs have the same index"""

    input_spec = SeparateTSInputSpec
    output_spec = SeparateTSOutputSpec

    def _run_interface(self, runtime):
            
        #import os
        #import numpy as np
        #import nibabel as nib
        
        #from neuropype_graph.utils_plot import plot_signals
        
        all_ts_file = self.inputs.all_ts_file
        
        path,fname_ts,ext = split_f(all_ts_file)
        
        ### loading ts shape = (trigs, electrods, time points)
        
        all_ts = np.load(all_ts_file)
        
        print "all_ts: " 
        print all_ts.shape
        
        separated_ts_files = []
        
        for i in range(all_ts.shape[0]):
            
            
            sep_ts_file = os.path.abspath(fname_ts + '_trig_' + str(i) + '.npy')
            
            np.save(sep_ts_file,all_ts[i,:,:])
            
            separated_ts_files.append(sep_ts_file)
        
        self.separated_ts_files = separated_ts_files
        
        return runtime
        
        #return mean_masked_ts_file,subj_coord_rois_file
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["separated_ts_files"] = self.separated_ts_files
        
        return outputs

################################################################################# RegressCovar ######################################################################################################################
 

#from neuropype_graph.utils_cor import regress_movement_wm_csf_parameters
from neuropype_graph.utils_cor import regress_parameters,regress_filter_normalize_parameters

class RegressCovarInputSpec(BaseInterfaceInputSpec):
    masked_ts_file = File(exists=True, desc='time series in npy format', mandatory=True)
    
    rp_file = File(exists=True, desc='Movement parameters', mandatory=True)
    
    mean_wm_ts_file = File(exists=True, desc='White matter signal', mandatory=False)
    
    mean_csf_ts_file = File(exists=True, desc='Cerebro-spinal fluid (ventricules) signal', mandatory=False)
    
    filtered = traits.Bool(True, usedefault = True , desc = "Is the signal filtered after regression?")
    
    normalized = traits.Bool(True, usedefault = True , desc = "Is the signal normalized after regression?")
    
class RegressCovarOutputSpec(TraitedSpec):
    
    resid_ts_file = File(exists=True, desc="residuals of time series after regression of all paramters")
    
class RegressCovar(BaseInterface):
    """
    Regress parameters of non-interest (i.e. movement parameters, white matter, csf) from signal
    Optionnally filter and normalize (z-score) the residuals
    """
    input_spec = RegressCovarInputSpec
    output_spec = RegressCovarOutputSpec

    def _run_interface(self, runtime):
                    
                    
        print "in regress_covariates"
        
        masked_ts_file = self.inputs.masked_ts_file
        rp_file = self.inputs.rp_file
        filtered = self.inputs.filtered
        normalized = self.inputs.normalized
        
        print "load masked_ts_file"
        
        data_mask_matrix = np.loadtxt(masked_ts_file)
        
        print data_mask_matrix.shape

        print "load rp parameters"
        
        print rp_file
        
        rp = np.genfromtxt(rp_file)
        #rp = np.loadtxt(rp_file,dtype = np.float)
        
        print rp.shape
        
        if isdefined(self.inputs.mean_csf_ts_file):
            
            mean_csf_ts_file = self.inputs.mean_csf_ts_file
            
            print "load mean_csf_ts_file" + str(mean_csf_ts_file)
            
            mean_csf_ts = np.loadtxt(mean_csf_ts_file)
            
            print mean_csf_ts.shape
            
            #rp = np.concatenate((rp,mean_csf_ts),axis = 1)
            rp = np.concatenate((rp,mean_csf_ts.reshape(mean_csf_ts.shape[0],1)),axis = 1)
            
            print rp.shape
            
            
        if isdefined(self.inputs.mean_wm_ts_file):
            
            mean_wm_ts_file = self.inputs.mean_wm_ts_file
            
            print "load mean_wm_ts_file"
            
            mean_wm_ts = np.loadtxt(mean_wm_ts_file)
            
            
            #rp = np.concatenate((rp,mean_csf_ts),axis = 1)
            rp = np.concatenate((rp,mean_wm_ts.reshape(mean_wm_ts.shape[0],1)),axis = 1)
            
            print rp.shape
            
        
        if filtered == True and normalized == True:
            
            ### regression movement parameters and computing z-score on the residuals
            #resid_data_matrix = regress_movement_wm_csf_parameters(data_mask_matrix,rp,mean_wm_ts,mean_csf_ts)
            resid_data_matrix,resid_filt_data_matrix,z_score_data_matrix = regress_filter_normalize_parameters(data_mask_matrix,rp)
            
            print resid_data_matrix.shape

            print "saving resid_ts"
            
            resid_ts_file = os.path.abspath('resid_ts.npy')
            np.save(resid_ts_file,z_score_data_matrix )

            print "plotting resid_ts"
            
            plot_resid_ts_file = os.path.abspath('resid_ts.eps')
            
            plot_sep_signals(plot_resid_ts_file,z_score_data_matrix)
            
            
            print "plotting diff filtered and non filtered data"
            
            plot_diff_filt_ts_file = os.path.abspath('diff_filt_ts.eps')
            
            plot_signals(plot_diff_filt_ts_file,np.array(resid_filt_data_matrix - resid_data_matrix,dtype = 'float'))
            
        elif filtered == False and normalized == False:
            
            print "Using only regression"
        
            ### regression movement parameters and computing z-score on the residuals
            #resid_data_matrix = regress_movement_wm_csf_parameters(data_mask_matrix,rp,mean_wm_ts,mean_csf_ts)
            resid_data_matrix = regress_parameters(data_mask_matrix,rp)
            
            print resid_data_matrix.shape

            print "saving resid_ts"
            
            resid_ts_file = os.path.abspath('resid_ts.npy')
            np.save(resid_ts_file,resid_data_matrix )

            
            resid_ts_txt_file = os.path.abspath('resid_ts.txt')
            np.savetxt(resid_ts_txt_file,resid_data_matrix,fmt = '%0.3f')

            
            print "plotting resid_ts"
            
            plot_resid_ts_file = os.path.abspath('resid_ts.eps')
            
            plot_sep_signals(plot_resid_ts_file,resid_data_matrix)
            
            
            #print "plotting diff filtered and non filtered data"
            
            #plot_diff_filt_ts_file = os.path.abspath('diff_filt_ts.eps')
            
            #plot_signals(plot_diff_filt_ts_file,np.array(resid_filt_data_matrix - resid_data_matrix,dtype = 'float'))
            
        else:
            
            print "Warning, not implemented (RegressCovar)"
            
        return runtime
        
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["resid_ts_file"] = os.path.abspath('resid_ts.npy')
    
        return outputs

        ################################################################################# FindSPMRegressor ######################################################################################################################
 

class FindSPMRegressorInputSpec(BaseInterfaceInputSpec):
    
    spm_mat_file = File(exists=True, desc='SPM design matrix after generate model', mandatory=True)
    
    regressor_name = traits.String(exists=True, desc='Name of the regressor in SPM design matrix to be looked after', mandatory=True)
    
    run_index = traits.Int(1, usedefault = True , desc = "Run (session) index, default is one in SPM")
     
    only_positive_values = traits.Bool(True, usedefault = True , desc = "Return only positive values of the regressor (negative values are set to 0)")
    
    concatenate_runs = traits.Int(1, usedefault = True , desc = "If concatenate runs, how many runs there is (needed to return the part of the regressors that is active for the session only)")
    
class FindSPMRegressorOutputSpec(TraitedSpec):
    
    regressor_file = File(exists=True, desc="txt file containing the regressor")
    
class FindSPMRegressor(BaseInterface):
    """
    Regress parameters of non-interest (i.e. movement parameters, white matter, csf) from signal
    Optionnally filter and normalize (z-score) the residuals
    """
    input_spec = FindSPMRegressorInputSpec
    output_spec = FindSPMRegressorOutputSpec

    def _run_interface(self, runtime):
                   
                   
                    
        import scipy.io
        import numpy as np
        import os

        #print spm_mat_file
        
        
        spm_mat_file = self.inputs.spm_mat_file
        regressor_name = self.inputs.regressor_name
        run_index = self.inputs.run_index
        only_positive_values = self.inputs.only_positive_values
        concatenate_runs = self.inputs.concatenate_runs
        
        
        print spm_mat_file
        
        ##Reading spm.mat for regressors extraction:
        d = scipy.io.loadmat(spm_mat_file)
        
        #print d
        
        
        ##Choosing the column according to the regressor name
        #_,col = np.where(d['SPM']['xX'][0][0]['name'][0][0] == u'Sn(1) ' + regressor_name)
        
        cond_name = u'Sn(' + str(run_index) + ') ' + regressor_name + '*bf(1)'
        
        print cond_name
        
        _,col = np.where(d['SPM']['xX'][0][0]['name'][0][0] == cond_name)
        
        print col
        
        ## reformating matrix (1,len) in vector (len)
        regressor_vect = d['SPM']['xX'][0][0]['X'][0][0][:,col].reshape(-1)
        

        print regressor_vect
        
        assert np.sum(regressor_vect) != 0, "Error, empty regressor {}".format(cond_name)
        
        if only_positive_values == True:
            
            regressor_vect[regressor_vect < 0] = 0
        
        if concatenate_runs != None:
            
            print run_index,concatenate_runs
            print regressor_vect.shape[0]
            
            nb_samples = regressor_vect.shape[0]/concatenate_runs
            
            print nb_samples
            
            begin_interval = (run_index-1)*nb_samples
            end_interval = run_index*nb_samples
            
            if 0 <= begin_interval and end_interval <= regressor_vect.shape[0]:
            
                print begin_interval,end_interval
            
                regressor_vect = regressor_vect[begin_interval:end_interval]
            
            else:
                
                print "Warning, error with interval [%d,%d]"%(begin_interval,end_interval)
        
            print regressor_vect.shape
            
        print "Saving extract_cond"
        regressor_file = os.path.abspath('extract_cond.txt')

        np.savetxt(regressor_file,regressor_vect)

        return runtime
        
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["regressor_file"] = os.path.abspath('extract_cond.txt')
        
        return outputs

        
        ################################################################################# MergeRuns ######################################################################################################################
 

class MergeRunsInputSpec(BaseInterfaceInputSpec):
    
    ts_files = traits.List(File(exists=True), desc='Numpy files with time series from different runs (sessions)',mandatory=True)

    regressor_files = traits.List(File(exists=True), desc='Txt files with regressors from different runs (sessions)', mandatory=True)
    
    coord_rois_files = traits.List(File(exists=True), desc='Txt files with coords from different runs (sessions)', mandatory=True)
    
class MergeRunsOutputSpec(TraitedSpec):
    
    ts_all_runs_file = File(exists=True, desc="npy file containing the merge ts")
    
    regressor_all_runs_file = File(exists=True, desc="txt file containing the merged regressors")
    
    coord_rois_all_runs_file = File(exists=True, desc="txt file containing the merged coords")
    
class MergeRuns(BaseInterface):
    """
    Merge time series,regressor files and coord files
    Could be done with different cases
    """
    input_spec = MergeRunsInputSpec
    output_spec = MergeRunsOutputSpec

    def _run_interface(self, runtime):
                   
        print 'in merge_runs'
        
        ts_files = self.inputs.ts_files
        regressor_files = self.inputs.regressor_files
        coord_rois_files = self.inputs.coord_rois_files
        
        
        if len(ts_files) != len(regressor_files):
            
            print "Warning, time series and regressors have different length (!= number of runs)"
            return 0
            
        if len(ts_files) != len(coord_rois_files):
            
            print "Warning, time series and number of coordinates have different length (!= number of runs)"
            return 0
            
        ### concatenate time series
        for i,ts_file in enumerate(ts_files):
            
            data_matrix = np.load(ts_file)
            
            print data_matrix.shape
            
            ## loading ROI coordinates
            coord_rois = np.loadtxt(coord_rois_files[i])
            
            print coord_rois.shape
            
            if i == 0:
                data_matrix_all_runs = np.empty((data_matrix.shape[0],0),dtype = data_matrix.dtype)
                
                coord_rois_all_runs = np.array(coord_rois,dtype = 'float')
                
                
            if coord_rois_all_runs.shape[0] != coord_rois.shape[0]:
                
                print "ROIs do not match for all different sessions "
                
                print os.getcwd()
                
                print "Warning, not implemented yet.... "
                
                ### pris de http://stackoverflow.com/questions/8317022/get-intersecting-rows-across-two-2d-numpy-arrays
                ### à tester....
                ### finir également la partie avec data_matrix_all_runs, en supprimant les colonnes qui ne sont pas communes à tous les runs...
                
                0/0
                
                A = coord_rois_all_runs
                B = coord_rois
                
                nrows, ncols = A.shape
                dtype={'names':['f{}'.format(i) for i in range(ncols)],
                    'formats':ncols * [A.dtype]}

                C = np.intersect1d(A.view(dtype), B.view(dtype))

                # This last bit is optional if you're okay with "C" being a structured array...
                C = C.view(A.dtype).reshape(-1, ncols)

                coord_rois_all_runs = C

            data_matrix_all_runs = np.concatenate((data_matrix_all_runs,data_matrix),axis = 1)
            
            print data_matrix_all_runs.shape
            
        ### save times series for all runs
        ts_all_runs_file = os.path.abspath('ts_all_runs.npy')
        
        np.save(ts_all_runs_file,data_matrix_all_runs)
        
        ### save coords in common for all runs
        coord_rois_all_runs_file = os.path.abspath('coord_rois_all_runs.txt')
        
        np.savetxt(coord_rois_all_runs_file,coord_rois_all_runs, fmt = '%2.3f')
        
        ### compute regressor for all sessions together (need to sum)
        
        print "compute regressor for all sessions together (need to sum)"
        
        regressor_all_runs = np.empty(shape = (0), dtype = float)
        
        ### Sum regressors
        for i,regress_file in enumerate(regressor_files):
            
            regress_data_vector = np.loadtxt(regress_file)
            
            if regress_data_vector.shape[0] != 0:
                
                if regressor_all_runs.shape[0] == 0:
                    
                    regressor_all_runs = regress_data_vector
                else:
                    regressor_all_runs = regressor_all_runs + regress_data_vector
                
            print np.sum(regressor_all_runs != 0.0)
        
        regressor_all_runs_file = os.path.abspath('regressor_all_runs.txt')

        np.savetxt(regressor_all_runs_file,regressor_all_runs,fmt = '%0.3f')

        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        print outputs
        
        outputs["ts_all_runs_file"] = os.path.abspath('ts_all_runs.npy')
        
        outputs["coord_rois_all_runs_file"] = os.path.abspath('coord_rois_all_runs.txt')
        
        outputs["regressor_all_runs_file"] = os.path.abspath('regressor_all_runs.txt')

        return outputs

        ################################################################################# ComputeConfCorMat ######################################################################################################################
 
from neuropype_graph.utils_cor import return_conf_cor_mat

from neuropype_graph.utils_plot import plot_hist,plot_cormat
        
class ComputeConfCorMatInputSpec(BaseInterfaceInputSpec):
    
    ts_file = File(exists=True, desc='Numpy files with time series to be correlated',mandatory=True)

    weight_file = File(exists=True, desc='Weight of the correlation (normally, condition regressor file)', mandatory=False)
    
    conf_interval_prob = traits.Float(0.05, usedefault = True, desc='Confidence interval', mandatory=True)
    
    plot_mat = traits.Bool(True, usedefault = True, desc='Confidence interval', mandatory=False)
    
    labels_file = File(exists=True, desc='Name of the nodes (used only if plot = true)', mandatory=False)
    
class ComputeConfCorMatOutputSpec(TraitedSpec):
    
    cor_mat_file = File(exists=True, desc="npy file containing the R values of correlation")
    
    Z_cor_mat_file = File(exists=True, desc="npy file containing the Z-values (after Fisher's R-to-Z trasformation) of correlation")
    
    conf_cor_mat_file = File(exists=True, desc="npy file containing the confidence interval around R values")
    
    #Z_conf_cor_mat_file = File(exists=True, desc="npy file containing the Z-values (after Fisher's R-to-Z trasformation) of correlation")
    
class ComputeConfCorMat(BaseInterface):
    """
    Compute correlation between time series, with a given confidence interval. If weight_file is specified, used for weighted correlation
    """
    input_spec = ComputeConfCorMatInputSpec
    output_spec = ComputeConfCorMatOutputSpec

    def _run_interface(self, runtime):
                   
                  
        print 'in compute_conf_correlation_matrix'
        
        ts_file = self.inputs.ts_file
        weight_file = self.inputs.weight_file
        conf_interval_prob = self.inputs.conf_interval_prob
            
        plot_mat = self.inputs.plot_mat
        labels_file = self.inputs.labels_file
        
        print 'load resid data'
        
        path, fname, ext = split_f(ts_file)
        
        data_matrix = np.load(ts_file)
        
        print data_matrix.shape
        
        print np.transpose(data_matrix).shape
        
        if isdefined(weight_file):
        
            print 'load weight_vect'
        
            weight_vect = np.loadtxt(weight_file)
            
            print weight_vect.shape
        
        else:
            weight_vect = np.ones(shape = (data_matrix.shape[1]))
        
        print "compute return_Z_cor_mat"
        
        cor_mat,Z_cor_mat,conf_cor_mat,Z_conf_cor_mat = return_conf_cor_mat(np.transpose(data_matrix),weight_vect,conf_interval_prob)
        
        #print cor_mat.shape
        
        #cor_mat = cor_mat + np.transpose(cor_mat)
        
        #sum_nan = np.sum(np.array(np.isnan(cor_mat), dtype = int),axis = 0)
        #sort_order =  np.argsort(sum_nan)
        
        #print sum_nan[sort_order]
        
        #print cor_mat[sort_order[-2],:]
        
        #a, = np.where(sum_nan == 1283)
        
        #print a
        
        #print a.shape
        #non_nan, = np.where(np.isnan(np.sum(cor_mat, axis = 0)))
        
        #print non_nan
        
        #nan, = np.where(np.sum(np.isnan(cor_mat), axis = 0) != 0)
        
        ##print nan
        
        #print "saving non_nan indexes as npy"
        
        #non_nan_file = os.path.abspath('non_nan_' + fname + '.npy')
        
        #np.save(non_nan_file,non_nan)
        
        
        #tmp = cor_mat[:,non_nan]
        
        #non_nan_cor_mat = tmp[non_nan,:]
        
        #new_non_nan = np.where(np.sum(np.isnan(non_nan_cor_mat), axis = 0) != 0)
        
        #print new_non_nan
        
        #0/0
        
        #non_nan_conf_cor_mat = conf_cor_mat[:,non_nan][non_nan,:]
        
        #non_nan_Z_cor_mat = Z_cor_mat[:,non_nan][non_nan,:]
        
        print "saving cor_mat as npy"
        
        cor_mat_file = os.path.abspath('cor_mat_' + fname + '.npy')
        
        #np.save(cor_mat_file,non_nan_cor_mat)
        np.save(cor_mat_file,cor_mat)
        
        print "saving conf_cor_mat as npy"
        
        conf_cor_mat_file = os.path.abspath('conf_cor_mat_' + fname + '.npy')
        
        #np.save(conf_cor_mat_file,non_nan_conf_cor_mat)
        np.save(conf_cor_mat_file,conf_cor_mat)
        
        print "saving Z_cor_mat as npy"
        
        Z_cor_mat_file = os.path.abspath('Z_cor_mat_' + fname + '.npy')
        
        #np.save(Z_cor_mat_file,non_nan_Z_cor_mat)
        np.save(Z_cor_mat_file,Z_cor_mat)
        
        
        #print "saving Z_conf_cor_mat as npy"
        
        #Z_conf_cor_mat_file = os.path.abspath('Z_conf_cor_mat_' + fname + '.npy')
        
        #np.save(Z_conf_cor_mat_file,Z_conf_cor_mat)
        
        
        if plot_mat == True:
            
            if isdefined(labels_file):
                    
                print 'extracting node labels'
                    
                labels = [line.strip() for line in open(labels_file)]
                print labels
                
            else:
                labels = []
            
            ############# cor_mat
            
            #### heatmap 
            
            print 'plotting cor_mat heatmap'
            
            plot_heatmap_cor_mat_file =  os.path.abspath('heatmap_cor_mat_' + fname + '.eps')
            
            plot_cormat(plot_heatmap_cor_mat_file,cor_mat,list_labels = labels)
            
            #### histogram 
            
            print 'plotting cor_mat histogram'
            
            plot_hist_cor_mat_file = os.path.abspath('hist_cor_mat_' + fname + '.eps')
            
            plot_hist(plot_hist_cor_mat_file,cor_mat,nb_bins = 100)
            
            ############ Z_cor_mat
            
            Z_cor_mat = np.load(Z_cor_mat_file)
            
            #### heatmap 
            
            print 'plotting Z_cor_mat heatmap'
            
            plot_heatmap_Z_cor_mat_file =  os.path.abspath('heatmap_Z_cor_mat_' + fname + '.eps')
            
            plot_cormat(plot_heatmap_Z_cor_mat_file,Z_cor_mat,list_labels = labels)
            
            #### histogram 
            
            print 'plotting Z_cor_mat histogram'
            
            plot_hist_Z_cor_mat_file = os.path.abspath('hist_Z_cor_mat_' + fname + '.eps')
            
            plot_hist(plot_hist_Z_cor_mat_file,Z_cor_mat,nb_bins = 100)
            
            ############ conf_cor_mat
            
            #### heatmap 
            
            print 'plotting conf_cor_mat heatmap'
            
            plot_heatmap_conf_cor_mat_file =  os.path.abspath('heatmap_conf_cor_mat_' + fname + '.eps')
            
            plot_cormat(plot_heatmap_conf_cor_mat_file,conf_cor_mat,list_labels = labels)
            
            #### histogram 
            
            print 'plotting conf_cor_mat histogram'
            
            plot_hist_conf_cor_mat_file = os.path.abspath('hist_conf_cor_mat_' + fname + '.eps')

            plot_hist(plot_hist_conf_cor_mat_file,conf_cor_mat,nb_bins = 100)
            
        
            ############# Z_conf_cor_mat
            
            #Z_conf_cor_mat = np.load(Z_conf_cor_mat_file)
            
            ##### heatmap 
            
            #print 'plotting Z_conf_cor_mat heatmap'
            
            #plot_heatmap_Z_conf_cor_mat_file =  os.path.abspath('heatmap_Z_conf_cor_mat_' + fname + '.eps')
            
            #plot_cormat(plot_heatmap_Z_conf_cor_mat_file,Z_conf_cor_mat,list_labels = labels)
            
            ##### histogram 
            
            #print 'plotting Z_conf_cor_mat histogram'
            
            #plot_hist_Z_conf_cor_mat_file = os.path.abspath('hist_Z_conf_cor_mat_' + fname + '.eps')
            
            #plot_hist(plot_hist_Z_conf_cor_mat_file,Z_conf_cor_mat,nb_bins = 100)
            
        
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        path, fname, ext = split_f(self.inputs.ts_file)
        
        outputs["cor_mat_file"] = os.path.abspath('cor_mat_' + fname + '.npy')
        
        outputs["conf_cor_mat_file"] = os.path.abspath('conf_cor_mat_' + fname + '.npy')
        
        outputs["Z_cor_mat_file"] = os.path.abspath('Z_cor_mat_' + fname + '.npy')
        
        #outputs["Z_conf_cor_mat_file"] = os.path.abspath('Z_conf_cor_mat_' + fname + '.npy')
        
        print outputs
        
        return outputs

        
        ################################################################################# SelectNonNAN ######################################################################################################################
 
class SelectNonNANInputSpec(BaseInterfaceInputSpec):
    
    sess_ts_files = traits.List(File(exists=True), desc='Numpy files with time series to be correlated',mandatory=True)

    sess_labels_files = traits.List(File(exists=True), desc='Name of the nodes (used only if plot = true)', mandatory=False)
    
class SelectNonNANOutputSpec(TraitedSpec):
    
    select_ts_files = traits.List(File(exists=True), desc='Numpy files with selected time series ',mandatory=True)

    select_labels_file = File(exists=True, desc="npy file containing the Z-values (after Fisher's R-to-Z trasformation) of correlation")
    
class SelectNonNAN(BaseInterface):
    
    """
    Select time series based on NaN
    TODO: finally no used
    """
    
    input_spec = SelectNonNANInputSpec
    output_spec = SelectNonNANOutputSpec

    def _run_interface(self, runtime):
                   
                  
        print 'in compute_conf_correlation_matrix'
        
        sess_ts_files = self.inputs.sess_ts_files
        labels_files = self.inputs.sess_labels_files
        
        
        if len(sess_ts_files) == 0:
            
            print "Warning, could not find sess_ts_files"
            
            return runtime
            
            
        path, fname, ext = split_f(sess_ts_files[0])
        
        list_sess_ts = []
        
        for ts_file in sess_ts_files:
            
            print 'load data'
            
            data_matrix = np.load(ts_file)
            
            print data_matrix.shape
            
            list_sess_ts.append(data_matrix)
            
        subj_ts = np.concatenate(tuple(list_sess_ts),axis = 0)
        
        print subj_ts.shape
        
        print np.sum(np.isnan(subj_ts) == True,axis = (1,2))
        print np.sum(np.isnan(subj_ts) == True,axis = (0,2))
        print np.sum(np.isnan(subj_ts) == True,axis = (0,1))
        
        good_trigs = np.sum(np.isnan(subj_ts) == True,axis = (1,2)) == 0
        
        select_subj_ts = subj_ts[good_trigs,:,:]
       
        print select_subj_ts.shape
        
        self.select_ts_files = []
        
        for i_trig in range(select_subj_ts.shape[0]):
        
            select_ts_file = os.path.abspath('select_' + fname + '_' + str(i_trig) + '.npy')
            
            np.save(select_ts_file,select_subj_ts[i_trig,:,:])
            
            self.select_ts_files.append(select_ts_file)
            
            
        ### check if all labels_files are identical
        
        if len(labels_files) == 0:
            
            print "Warning, could not find sess_ts_files"
            
            return runtime
          
        select_labels_file = labels_files[0]
        
        select_labels = np.array(np.loadtxt(select_labels_file),dtype = 'str')
        
        print select_labels
        0/0
        labels_files
        
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["select_ts_files"] = self.select_ts_files
        
        print outputs
        
        return outputs

############################################# used in run_mean_correl ###############################################################################

##################################################### PrepareMeanCorrel ###############################################################################

from neuropype_graph.utils_cor import return_corres_correl_mat

class PrepareMeanCorrelInputSpec(BaseInterfaceInputSpec):
    
    #print "in PrepareMeanCorrelInputSpec"
    
    gm_mask_coords_file = File(exists=True, desc='reference coordinates',mandatory=True)
                
    cor_mat_files = traits.List(File(exists=True), desc='Numpy files with correlation matrices gm_mask_coords',mandatory=True)

    coords_files = traits.List(File(exists=True), desc='Txt files with coordinates (corresponding to the space also described in ', mandatory=True)
    
    labels_file = File(exists=True, desc='reference labels',mandatory=False)
    
    plot_mat = traits.Bool(True, usedefault = True, mandatory = False)
    
class PrepareMeanCorrelOutputSpec(TraitedSpec):
    
    #print "in PrepareMeanCorrelOutputSpec"
    
    group_cor_mat_matrix_file = File(exists=True, desc="npy file containing all correlation matrices in 3D")
    
    sum_cor_mat_matrix_file = File(exists=True, desc="npy file containing the sum of all correlation matrices")
    
    sum_possible_edge_matrix_file = File(exists=True, desc="npy file containing the number of correlation matrices where both nodes where actually defined in the mask")
    
    avg_cor_mat_matrix_file = File(exists=True, desc="npy file containing the average of all correlation matrices")
    
class PrepareMeanCorrel(BaseInterface):
    
    import numpy as np
    import os

    #import nibabel as nib
    
    """
    Return average of correlation values within the same common space (defined in gm_mask_coords), only when the nodes are defined for a given values 
    """
    
    input_spec = PrepareMeanCorrelInputSpec
    output_spec = PrepareMeanCorrelOutputSpec

    def _run_interface(self, runtime):
                
            gm_mask_coords_file = self.inputs.gm_mask_coords_file
            cor_mat_files = self.inputs.cor_mat_files
            coords_files = self.inputs.coords_files
            labels_file  = self.inputs.labels_file
            plot_mat  = self.inputs.plot_mat
            
            print 'loading gm mask corres'
            
            gm_mask_coords = np.loadtxt(gm_mask_coords_file)
            
            print gm_mask_coords.shape
                
            #### read matrix from the first group
            #print Z_cor_mat_files
            
            sum_cor_mat_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = float)
            print sum_cor_mat_matrix.shape
            
            sum_possible_edge_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = int)
            print sum_possible_edge_matrix.shape
            
                    
            group_cor_mat_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0],len(cor_mat_files)),dtype = float)
            print group_cor_mat_matrix.shape
            
            if len(cor_mat_files) != len(coords_files):
                print "warning, length of cor_mat_files, coords_files are imcompatible {} {} {}".format(len(cor_mat_files),len(coords_files))
            
            for index_file in range(len(cor_mat_files)):
                
                print cor_mat_files[index_file]
                
                if os.path.exists(cor_mat_files[index_file]) and os.path.exists(coords_files[index_file]):
                
                    Z_cor_mat = np.load(cor_mat_files[index_file])
                    print Z_cor_mat.shape
                    
                    
                    coords = np.loadtxt(coords_files[index_file])
                    #print coords.shape
                    
                    
                    
                    corres_cor_mat,possible_edge_mat = return_corres_correl_mat(Z_cor_mat,coords,gm_mask_coords)
                    
                    
                    np.fill_diagonal(corres_cor_mat,0)
                    
                    np.fill_diagonal(possible_edge_mat,1)
                    
                    sum_cor_mat_matrix += corres_cor_mat
                    
                    sum_possible_edge_matrix += possible_edge_mat
                    
                    
                    group_cor_mat_matrix[:,:,index_file] = corres_cor_mat
                    
                    
                else:
                    print "Warning, one or more files between " + cor_mat_files[index_file] + ', ' + coords_files[index_file] + " do not exists"
                
                
            self.group_cor_mat_matrix_file= os.path.abspath('group_cor_mat_matrix.npy')
            
            np.save(self.group_cor_mat_matrix_file,group_cor_mat_matrix)
            
            
            print 'saving sum cor_mat matrix'
            
            self.sum_cor_mat_matrix_file = os.path.abspath('sum_cor_mat_matrix.npy')
            
            np.save(self.sum_cor_mat_matrix_file,sum_cor_mat_matrix)
            
            print 'saving sum_possible_edge matrix'
            
            self.sum_possible_edge_matrix_file = os.path.abspath('sum_possible_edge_matrix.npy')
            
            np.save(self.sum_possible_edge_matrix_file,sum_possible_edge_matrix)
            
            print 'saving avg_cor_mat_matrix'
            
            self.avg_cor_mat_matrix_file  = os.path.abspath('avg_cor_mat_matrix.npy')
            
            avg_cor_mat_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = float)
            
            
            if (np.sum(np.array(sum_possible_edge_matrix == 0)) != 0):
            
                    avg_cor_mat_matrix = np.divide(np.array(sum_cor_mat_matrix,dtype = float),np.array(sum_possible_edge_matrix,dtype = float))
                            
                    avg_cor_mat_matrix[np.isnan(avg_cor_mat_matrix)] = 0.0
                    
                    print np.amin(avg_cor_mat_matrix),np.amax(avg_cor_mat_matrix)
                    
                    
                    np.save(self.avg_cor_mat_matrix_file,avg_cor_mat_matrix)
            
            else:
                    print "!!!!!!!!!!!!!!!!!!!!!!Breaking!!!!!!!!!!!!!!!!!!!!!!!!, found 0 elements in sum_cor_mat_matrix"
                    return
                
                
            if plot_mat == True:
                
                if isdefined(labels_file):
                        
                    print 'extracting node labels'
                        
                    labels = [line.strip() for line in open(labels_file)]
                    print labels
                    
                else:
                    labels = []
                
                #### heatmap 
                
                print 'plotting Z_cor_mat heatmap'
                
                plot_heatmap_avg_cor_mat_file =  os.path.abspath('heatmap_avg_cor_mat.eps')
                
                plot_cormat(plot_heatmap_avg_cor_mat_file,avg_cor_mat_matrix,list_labels = labels)
                
                
            return runtime
                    

    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["group_cor_mat_matrix_file"] = self.group_cor_mat_matrix_file
        outputs["sum_cor_mat_matrix_file"] = self.sum_cor_mat_matrix_file
        outputs["sum_possible_edge_matrix_file"] = self.sum_possible_edge_matrix_file
        outputs["avg_cor_mat_matrix_file"] = self.avg_cor_mat_matrix_file
        
        #print outputs
        
        return outputs



##################################################### PreparePermutMeanCorrel ###############################################################################

class PreparePermutMeanCorrelInputSpec(BaseInterfaceInputSpec):
       
    cor_mat_files = traits.List(File(exists=True), desc='Numpy files with correlation matrices gm_mask_coords',mandatory=True)

    permut_group_sizes = traits.List(traits.Int, desc = 'How to split the groups after shuffling', mandatory = True)
    
    seed = traits.Int(0, usedefault = True, decs = 'Start of random process')
    
    #variance_adjst = traits.Bool(False, usedefault = True, desc = "is between-subject variance adjusted taken into account?")
    
    
class PreparePermutMeanCorrelOutputSpec(TraitedSpec):
    
    permut_mean_cormat_files = traits.List(File(exists=True),desc="npy files containing the average of permuted correlation matrices")
    
class PreparePermutMeanCorrel(BaseInterface):
    
    import numpy as np
    import os

    #import nibabel as nib
    
    """
    Return average of correlation values after shuffling orig datasets
    """
    
    input_spec = PreparePermutMeanCorrelInputSpec
    output_spec = PreparePermutMeanCorrelOutputSpec

    def _run_interface(self, runtime):
               
            print self.inputs.seed
            
            np.random.seed(self.inputs.seed)
            
            cormats = [np.load(cor_mat_file) for cor_mat_file in self.inputs.cor_mat_files]
            
            print cormats
            
            assert len(cormats) == sum(self.inputs.permut_group_sizes), "Warning, len(cormats) {0} != sum permut_group_sizes {1}".format(len(cormats),sum(self.inputs.permut_group_sizes))
            
            subj_indexes = np.arange(len(cormats))
            
            np.random.shuffle(subj_indexes)
            
            print subj_indexes
            
            subj_indexes_file = os.path.abspath("subj_indexes.txt")

            f = open(subj_indexes_file,"w+")
            
            #f.write(subj_indexes.tolist())
            
            np.savetxt(f,subj_indexes,fmt = "%d")
                  
                  
            min_index = 0
            
            cormats = np.array(cormats)
            
            print cormats.shape
            
            self.permut_mean_cormat_files = []
            
            for i,cur_nb_subj_by_gender in enumerate(self.inputs.permut_group_sizes):
                
                print cur_nb_subj_by_gender
                
                cur_range = np.arange(min_index,min_index+cur_nb_subj_by_gender)
                
                print cur_range
                
                rand_indexes = subj_indexes[cur_range]
                
                print rand_indexes
                
                #f.write(rand_indexes)
                
                np.savetxt(f,rand_indexes,fmt = "%d")
                
                permut_cormats = cormats[rand_indexes,:,:]
                
                permut_mean_cormat = np.mean(permut_cormats,axis = 0)
                
                print permut_mean_cormat.shape
                
                permut_mean_cormat_file = os.path.abspath("permut_mean_cormat_" + str(i) + ".npy")
                
                np.save(permut_mean_cormat_file,permut_mean_cormat)
                
                self.permut_mean_cormat_files.append(permut_mean_cormat_file)
                
                min_index+=cur_nb_subj_by_gender
            f.close()
            
            return runtime
                    

    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["permut_mean_cormat_files"] = self.permut_mean_cormat_files
        
        #print outputs
        
        return outputs

                    

#def prepare_signif_correlation_matrices(cor_mat_files,conf_cor_mat_files,coords_files,gm_mask_coords_file):


    #import numpy as np
    #import os

    ##import nibabel as nib
    
    ##from utils_cor import read_Pajek_corres_nodes,read_lol_file
    
    #print 'loading gm mask corres'
    
    #gm_mask_coords = np.loadtxt(gm_mask_coords_file)
    
    #print gm_mask_coords.shape
        
    ##### read matrix from the first group
    ##print Z_cor_mat_files
    
    #sum_signif_cor_mat = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = float)
    #print sum_signif_cor_mat.shape
    
    #sum_possible_edge_mat = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = int)
    #print sum_possible_edge_mat.shape
    
    
    
    #group_signif_cor_mat = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0],len(cor_mat_files)),dtype = float)
    #print group_signif_cor_mat.shape
    
    #if len(cor_mat_files) != len(coords_files):
        #print "warning, length of cor_mat_files, coords_files are imcompatible {} {} {}".format(len(cor_mat_files),len(coords_files))
    
    #for index_file in range(len(cor_mat_files)):
        
        #print cor_mat_files[index_file]
        
        #if os.path.exists(cor_mat_files[index_file]) and os.path.exists(coords_files[index_file]) and os.path.exists(conf_cor_mat_files[index_file]):
        
            #cor_mat = np.load(cor_mat_files[index_file])
            #print cor_mat.shape
            
            #conf_cor_mat = np.load(conf_cor_mat_files[index_file])
            #print conf_cor_mat.shape
            
            #signif_cor_mat = conf_cor_mat < np.abs(cor_mat)
            
            #print signif_cor_mat.shape
            
            #coords = np.loadtxt(coords_files[index_file])
            #print coords.shape
            
            #corres_signif_cor_mat,possible_edge_mat = return_corres_correl_mat(signif_cor_mat,coords,gm_mask_coords)
            
            #group_signif_cor_mat[:,:,index_file] = corres_signif_cor_mat
            
            #sum_signif_cor_mat += corres_signif_cor_mat
            
            #sum_possible_edge_mat += possible_edge_mat
            
        #else:
            #print "Warning, one or more files between " + cor_mat_files[index_file] + ', ' + coords_files[index_file] + ', ' + conf_cor_mat_files[index_file]+ " do not exists"
        
    #print 'computing sum_signif_cor_mat'
    
    #sum_signif_cor_mat = np.sum(group_signif_cor_mat,axis = 2)
    
    #print 'saving group_signif cor_mat matrix'
    
    #group_signif_cor_mat_file= os.path.abspath('group_signif_cor_mat.npy')
    
    #np.save(group_signif_cor_mat_file,group_signif_cor_mat)
    
    
    #print 'saving sum_signif cor_mat matrix'
    
    #sum_signif_cor_mat_file = os.path.abspath('sum_signif_cor_mat.npy')
    
    #np.save(sum_signif_cor_mat_file,sum_signif_cor_mat)
    
    
    #print 'saving sum_possible_edge matrix'
    
    #sum_possible_edge_mat_file = os.path.abspath('sum_possible_edge_mat.npy')
    
    #np.save(sum_possible_edge_mat_file,sum_possible_edge_mat)
    
    
    
    
    
    #print 'saving avg_cor_mat_matrix'
    
    #norm_signif_cor_mat_file  = os.path.abspath('norm_signif_cor_mat.npy')
    
    #norm_signif_cor_mat = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = int)
    
    #if (np.where(np.array(sum_possible_edge_mat == 0)) != 0):
    
            #norm_signif_cor_mat = np.divide(np.array(sum_signif_cor_mat,dtype = float),np.array(sum_possible_edge_mat,dtype = float))
            
            #np.save(norm_signif_cor_mat_file,norm_signif_cor_mat)
    
    #else:
            #print "!!!!!!!!!!!!!!!!!!!!!!Breaking!!!!!!!!!!!!!!!!!!!!!!!!, found 0 elements in norm_signif_cor_mat"
            #return
            
    
    #return group_signif_cor_mat_file,sum_signif_cor_mat_file,sum_possible_edge_mat_file,norm_signif_cor_mat_file
        
     
        
        
