# -*- coding: utf-8 -*-
"""
Support for computing ROIs mask by different means:
- peak activations
- MNI coordinates
- ROI files (in nifti format)
- from template (Harvard Oxford = HO)

The outputs of the functions will always be a labeled mask, with values starting from 1 (0 being the background image)

"""

import nipype.interfaces.spm as spm

from nipype.utils.filemanip import split_filename as split_f

from neuropype_graph.utils import check_np_dimension

import itertools as iter

from nipy.labs.viz import coord_transform

import numpy as np
import nibabel as nib
import glob
import os


from xml.dom import minidom
############################################################ from a list of MNI coords ############################################################################

### utils for merging different label and coord file before computing label masks if necessary
def merge_coord_and_label_files(ROI_coords_dir):
    """
    utils for merging different label and coord file before computing label masks if necessary
    should be rewritten to more more general with glob("Coord*.txt) and glob("Labels*.txt)...
    """

    import os
    import numpy as np
    
    list_coords = []
    
    list_labels = []
        
    for event in ['Odor','Recall']:
    
        print event
        
        event_coords_file = os.path.join(ROI_coords_dir,"Coord_Network"+event+".txt")
        
        print event_coords_file
        
        for line in open(event_coords_file):
        
            print line
            
            list_coord = map(int,line.strip().split('\t'))
            
            print list_coord
            
            list_coords.append(list_coord)
            
            
        event_labels_file = os.path.join(ROI_coords_dir,"Labels_Network"+event+".txt")
        
        for line in open(event_labels_file):
        
            print line
            
            list_labels.append(line.strip())
                    
        
    print list_coords
    
    print list_labels
        
    print len(list_coords)
    print len(list_labels)
    
    ###### saving merged file
    all_coords_file = os.path.join(ROI_coords_dir,"Coord_Network_Odor-Recall.txt")
    all_labels_file = os.path.join(ROI_coords_dir,"Labels_Network_Odor-Recall.txt")
    
    np.savetxt(all_coords_file,np.array(list_coords,dtype = 'int'),fmt = '%d')
    np.savetxt(all_labels_file,np.array(list_labels,dtype = 'string'),fmt = '%s')
    
    return all_coords_file,all_labels_file
    
def compute_labelled_mask_from_ROI_coords_files(ref_img_file,MNI_coords_file,neighbourhood = 1):
    """
    Compute labeled mask by specifying MNI coordinates and labels 'at hand'
    #"""
    
    ref_image = nib.load(ref_img_file)
    
    ref_image_data = ref_image.get_data()
    
    ref_image_data_shape = ref_image_data.shape
    
    print ref_image_data_shape
    
    ref_image_data_sform = ref_image.get_sform()
    
    print ref_image_data_sform
    
    ROI_MNI_coords_list = np.array(np.loadtxt(MNI_coords_file),dtype = 'int').tolist()
    
    print ROI_MNI_coords_list
    
    #ROI_labels = [lign.strip() for lign in open(labels_file)]
    
    #print labels
    
    print len(ROI_MNI_coords_list)
    #print len(ROI_labels)
    
    ### transform MNI coords to numpy coords
    ## transfo inverse de celle stockes dans le header
    mni_sform_inv = np.linalg.inv(ref_image_data_sform)
    
    ROI_coords = np.array([coord_transform(x, y, z, mni_sform_inv) for x,y,z in ROI_MNI_coords_list],dtype = "int64")
    
    ROI_coords_labelled_mask = np.zeros(shape = ref_image_data_shape,dtype = 'int64') - 1
    
    print ROI_coords_labelled_mask
    
    for i,ROI_coord in enumerate(ROI_coords):
    
        print ROI_coord
        
        for relative_coord in iter.product(range(-neighbourhood,neighbourhood+1), repeat=3):

            neigh_x,neigh_y,neigh_z = ROI_coord + relative_coord

            print neigh_x,neigh_y,neigh_z
            
            if check_np_dimension(ROI_coords_labelled_mask.shape,np.array([neigh_x,neigh_y,neigh_z],dtype = 'int64')):
            
                ROI_coords_labelled_mask[neigh_x,neigh_y,neigh_z] = i
            
           
    print np.unique(ROI_coords_labelled_mask)
    
    #path, fname, ext = '','',''
    path, fname, ext = split_f(MNI_coords_file)
    
    #ROI_coords_labelled_mask_file = os.path.join(path,"All_labelled_ROI2-neigh_"+str(neighbourhood)+".nii")
    ROI_coords_labelled_mask_file = os.path.join(path,"All_labelled_ROI-neigh_"+str(neighbourhood)+".nii")
    
    #ROI_coords_np_coords_file = os.path.join(path,"All_ROI_np_coords.txt")
    ROI_coords_np_coords_file = os.path.join(path,"All_ROI_np_coords.txt")
    
    ###save ROI_coords_labelled_mask
    nib.save(nib.Nifti1Image(ROI_coords_labelled_mask,ref_image.get_affine(),ref_image.get_header()),ROI_coords_labelled_mask_file)
    
    #### save np coords
    np.savetxt(ROI_coords_np_coords_file,np.array(ROI_coords,dtype = int),fmt = "%d")
    
    
    return ROI_coords_labelled_mask_file
    
############################################################ from a list of MNI coords (output one VOI binary mask nii image) ############################################################################
def compute_ROI_nii_from_ROI_coords_files(ref_img_file,MNI_coords_file,labels_file,neighbourhood = 1):
    """
    Export single file VOI binary nii image 
    #"""
    
    ref_image = nib.load(ref_img_file)
    
    ref_image_data = ref_image.get_data()
    
    ref_image_data_shape = ref_image_data.shape
    
    print ref_image_data_shape
    
    ref_image_data_sform = ref_image.get_sform()
    
    print ref_image_data_sform
    
    ROI_MNI_coords_list = np.array(np.loadtxt(MNI_coords_file),dtype = 'int').tolist()
    
    print ROI_MNI_coords_list
    
    ROI_labels = [lign.strip() for lign in open(labels_file)]
    
    print ROI_labels
    
    print len(ROI_MNI_coords_list)
    print len(ROI_labels)
        
    ### transform MNI coords to numpy coords
    ## transfo inverse de celle stockes dans le header
    mni_sform_inv = np.linalg.inv(ref_image_data_sform)
    
    ROI_coords = np.array([coord_transform(x, y, z, mni_sform_inv) for x,y,z in ROI_MNI_coords_list],dtype = "int64")
        
    for i,ROI_coord in enumerate(ROI_coords):
    
        ROI_coords_labelled_mask = np.zeros(shape = ref_image_data_shape,dtype = 'int64')
    
        print ROI_coord
        print ROI_labels[i]
        
        for relative_coord in iter.product(range(-neighbourhood,neighbourhood+1), repeat=3):

            neigh_x,neigh_y,neigh_z = ROI_coord + relative_coord

            print neigh_x,neigh_y,neigh_z
            
            if check_np_dimension(ROI_coords_labelled_mask.shape,np.array([neigh_x,neigh_y,neigh_z],dtype = 'int64')):
            
                ROI_coords_labelled_mask[neigh_x,neigh_y,neigh_z] = 1
            
            
        print ROI_coords_labelled_mask
    
        #path, fname, ext = '','',''
        path, fname, ext = split_f(MNI_coords_file)
        
        
        ROI_coords_labelled_mask_file = os.path.join(path,"ROI_" + ROI_labels[i] + "-neigh_" + str(neighbourhood) + "_2.nii")
            
        ###save ROI_coords_labelled_mask
        nib.save(nib.Nifti1Image(ROI_coords_labelled_mask,ref_image.get_affine(),ref_image.get_header()),ROI_coords_labelled_mask_file)
        
    return ROI_coords_labelled_mask_file

#### from a ROI directories, containing a list of VOI binary mask nii images ####
def compute_labelled_mask_from_anat_ROIs(ref_img_file,ROI_dir):
    
    """
    compute labelled_mask from a list of img files, presenting ROIs extracted from MRIcron in the nii or img format
    each ROI is represented by a different IMG file and should start by 'ROI_'. Resampling is done based on the shape of ref_img_file
    """
    
    ref_image = nib.load(ref_img_file)
    
    ref_image_data = ref_image.get_data()
    
    ref_image_data_shape = ref_image_data.shape
    
    print ROI_dir
    
    resliced_ROI_files =  glob.glob(os.path.join(ROI_dir,"rROI*.nii"))

    print resliced_ROI_files
    print len(resliced_ROI_files)
    
    ROI_files =  glob.glob(os.path.join(ROI_dir,"ROI*.nii"))
    
    print ROI_files
    print len(ROI_files)
    
    if len(resliced_ROI_files) != len(ROI_files) :
        
        for i,ROI_file in enumerate(ROI_files):
        
            ROI_image = nib.load(ROI_file)
        
            ROI_data = ROI_image.get_data()
            
            ROI_data_shape = ROI_data.shape
            
            print "Original ROI template %d shape:"%i
            
            print ROI_data.shape
            
            reslice_ROI = spm.Reslice()
            reslice_ROI.inputs.in_file = ROI_file
            reslice_ROI.inputs.space_defining = ref_img_file

            resliced_ROI_file =  reslice_ROI.run().outputs.out_file

        resliced_ROI_files = glob.glob(os.path.join(ROI_dir,"rROI*.nii"))
        
    resliced_ROI_files.sort()
    
    print resliced_ROI_files
    print len(resliced_ROI_files)
    
    labels = []
    
    labelled_mask_data = np.zeros(shape = ref_image_data.shape, dtype = 'int')
    
    print labelled_mask_data.shape
    
    for i,resliced_ROI_file in enumerate(resliced_ROI_files):
    
        print i 
        
        path,fname,ext = split_f(resliced_ROI_file)
    
        labels.append(fname)
    
        resliced_ROI_img = nib.load(resliced_ROI_file)

        resliced_ROI_data = resliced_ROI_img.get_data()

        print resliced_ROI_data.shape
        
        print np.sum(resliced_ROI_data != 0)
        
        labelled_mask_data[resliced_ROI_data != 0] = i+1
        
        print np.unique(labelled_mask_data)
    
    print np.unique(labelled_mask_data).shape
    
    print len(labels)
    
    ### save labeled_mask
    labelled_mask_data_file = os.path.join(ROI_dir,"all_ROIs_labelled_mask.nii")
    
    nib.save(nib.Nifti1Image(labelled_mask_data,ref_image.get_affine(),ref_image.get_header()),labelled_mask_data_file)
    
    ### save labels
    labels_list_file = os.path.join(ROI_dir,"labels_all_ROIs.txt")
    
    np.savetxt(labels_list_file,np.array(labels,dtype = 'string'),fmt = "%s")
    
    return labelled_mask_data_file,labels_list_file
    #nib.load(ref_img_file)
    
def compute_MNI_coords_from_indexed_template(indexed_template_file,ROI_dir):
    
    """
    compute MNI coords from an indexed template
    """
    
    path, base,ext = split_f(indexed_template_file)
    
    print base
    
    if len(base.split("-")) > 1:
        base_name = base.split("-")[1]
    else:
        base_name = base
          
    ref_image = nib.load(indexed_template_file)
    
    ref_image_data = ref_image.get_data()
    
    print ref_image_data.shape
    
    ref_image_affine = ref_image.affine
    
    print ref_image_affine
    
    ROI_coords = []
    
    ROI_MNI_coords = []
    
    for index in np.unique(ref_image_data):
        
        i,j,k = np.where(ref_image_data == index)
        
        mean_coord_ijk = np.mean(np.array((i,j,k)), axis = 1)
        
        print mean_coord_ijk
        
        ROI_coords.append(mean_coord_ijk)
        
        MNI_coord = np.dot(ref_image_affine,np.append(mean_coord_ijk,1))
        
        print MNI_coord
        
        ROI_MNI_coords.append(MNI_coord[:3])
        
    ROI_coords = np.array(ROI_coords,dtype = float)
    
    print ROI_coords
    
    ROI_MNI_coords = np.array(ROI_MNI_coords,dtype = float)
    
    print ROI_MNI_coords
    
    ROI_coords_file = os.path.join(ROI_dir,"ROI_coords-" + base_name + ".txt")
    
    np.savetxt(ROI_coords_file,ROI_coords, fmt = "%.3f %.3f %.3f")
    
    
    ROI_MNI_coords_file = os.path.join(ROI_dir,"ROI_MNI_coords-" + base_name + ".txt")
    
    np.savetxt(ROI_MNI_coords_file,ROI_MNI_coords, fmt = "%.3f %.3f %.3f")
    
    
    return ROI_coords_file,ROI_MNI_coords_file
    #nib.load(ref_img_file)


################################# preparing HO template by recombining sub and cortical mask + reslicing to image format ##################################

def compute_recombined_HO_template(img_file,ROI_dir,HO_dir = "/usr/share/fsl/data/atlases/"):


    if not os.path.exists(ROI_dir):
        
        os.makedirs(ROI_dir)
        
    #img = nib.load(img_file)
    
    #img_header = img.get_header()
    
    #img_affine = img.get_affine()
    
    #img_shape = img.get_data().shape

    #print img_shape
    ### cortical
    HO_cortl_img_file = os.path.join(HO_dir,"HarvardOxford/HarvardOxford-cortl-maxprob-thr25-2mm.nii.gz")
    
    HO_cortl_img = nib.load(HO_cortl_img_file)
    
    HO_cortl_data = HO_cortl_img.get_data()
    
    print HO_cortl_data.shape
    
    print np.min(HO_cortl_data),np.max(HO_cortl_data)
    
    ### subcortical
    HO_sub_img_file = os.path.join(HO_dir,"HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz")
    
    HO_sub_img = nib.load(HO_sub_img_file)
    
    HO_sub_img_header = HO_sub_img.get_header()
    
    HO_sub_img_affine = HO_sub_img.get_affine()
     
    HO_sub_data = HO_sub_img.get_data()
    
    print HO_sub_data.shape
    
    print np.min(HO_sub_data),np.max(HO_sub_data)
    
    ################# White matter mask
    white_matter_HO_img_file = os.path.join(ROI_dir, "Harvard-Oxford-white_matter.nii")
    
    if not os.path.isfile(white_matter_HO_img_file):
        
        #### extracting mask for white matter :
        white_matter_HO_data = np.zeros(shape = HO_sub_data.shape,dtype = 'int')
        
        ###left white matter
        white_matter_HO_data[HO_sub_data == 1] = 1
        
        ### right white matter
        white_matter_HO_data[HO_sub_data == 12] = 1
        
        #### saving white matter mask
        white_matter_HO_data = np.array(white_matter_HO_data,dtype = 'int')
        
        nib.save(nib.Nifti1Image(dataobj = white_matter_HO_data,header = HO_sub_img_header,affine = HO_sub_img_affine),white_matter_HO_img_file)
        
    else:
        
        white_matter_HO_data = nib.load(white_matter_HO_img_file).get_data()
        
    ###reslicing to ref image
    resliced_white_matter_HO_img_file = os.path.join(ROI_dir, "rHarvard-Oxford-white_matter.nii")
    
    if not os.path.isfile(resliced_white_matter_HO_img_file):
    
        reslice_white_matter_HO = spm.Reslice()
        reslice_white_matter_HO.inputs.in_file = white_matter_HO_img_file
        reslice_white_matter_HO.inputs.space_defining = img_file
        reslice_white_matter_HO.inputs.out_file = resliced_white_matter_HO_img_file
        
        reslice_white_matter_HO.run()
        
    ################# grey matter mask
    grey_matter_HO_img_file = os.path.join(ROI_dir, "Harvard-Oxford-grey_matter.nii")
    
    if not os.path.isfile(grey_matter_HO_img_file):
        
        #### extracting mask for grey matter :
        grey_matter_HO_data = np.zeros(shape = HO_sub_data.shape,dtype = 'int')
        
        ###left grey matter
        grey_matter_HO_data[HO_sub_data == 2] = 1
        
        ### right grey matter
        grey_matter_HO_data[HO_sub_data == 13] = 1
        
        #### saving grey matter mask
        grey_matter_HO_data = np.array(grey_matter_HO_data,dtype = 'int')
        nib.save(nib.Nifti1Image(dataobj = grey_matter_HO_data,header = HO_sub_img_header,affine = HO_sub_img_affine),grey_matter_HO_img_file)
        
    else:
        grey_matter_HO_data = nib.load(grey_matter_HO_img_file).get_data()
        
    ### reslicing to ref image
    resliced_grey_matter_HO_img_file = os.path.join(ROI_dir, "rHarvard-Oxford-grey_matter.nii")
    
    if not os.path.isfile(resliced_grey_matter_HO_img_file):
    
        reslice_grey_matter_HO = spm.Reslice()
        reslice_grey_matter_HO.inputs.in_file = grey_matter_HO_img_file
        reslice_grey_matter_HO.inputs.space_defining = img_file
        reslice_grey_matter_HO.inputs.out_file = resliced_grey_matter_HO_img_file
        
        reslice_grey_matter_HO.run()
        
    
        
    ######## Ventricule (+ apprently outside the brain) mask
    ventricule_HO_img_file = os.path.join(ROI_dir, "Harvard-Oxford-ventricule.nii")
    
    if not os.path.isfile(ventricule_HO_img_file):
        
        #### extracting mask for ventricules:
        ventricule_HO_data = np.zeros(shape = HO_sub_data.shape,dtype = 'int')
        
        ### left ventricule
        ventricule_HO_data[HO_sub_data == 3] = 1
        
        ### right ventricule
        ventricule_HO_data[HO_sub_data == 14] = 1
        
        #### saving white matter mask
        ventricule_HO_data = np.array(ventricule_HO_data,dtype = 'int')
        nib.save(nib.Nifti1Image(dataobj = ventricule_HO_data,header = HO_sub_img_header,affine = HO_sub_img_affine),ventricule_HO_img_file)
          
    else:
        ventricule_HO_data = nib.load(ventricule_HO_img_file).get_data()
         
    ### reslicing to ref image
    resliced_ventricule_HO_img_file = os.path.join(ROI_dir, "rHarvard-Oxford-ventricule.nii")
    
    if not os.path.isfile(resliced_ventricule_HO_img_file):
    
        reslice_ventricule_HO = spm.Reslice()
        reslice_ventricule_HO.inputs.in_file = ventricule_HO_img_file
        reslice_ventricule_HO.inputs.space_defining = img_file
        reslice_ventricule_HO.inputs.out_file = resliced_ventricule_HO_img_file
        
        reslice_ventricule_HO.run()
        
    
    ################################ merging data from cortl and sub 
    #useful_sub_indexes = [1] + range(3,11) + [12] + range(14,21)
    
    ### sans BrainStem
    useful_sub_indexes = [1] + range(3,7) + range(8,11) + [12] + range(14,21)
    print useful_sub_indexes
        
    useful_cortl_indexes = np.unique(HO_cortl_data)[:-1]
    
    print useful_cortl_indexes
    
    ####### concatenate areas from cortical and subcortical masks
    full_HO_img_file = os.path.join(ROI_dir, "Harvard-Oxford-cortl-sub-recombined-111regions.nii")
    
    if not os.path.isfile(full_HO_img_file):
        
        ### recombining indexes
        full_HO_data = np.zeros(shape = HO_cortl_data.shape,dtype = 'int')
        
        #print full_HO_data
        
        
        new_index = 1
        
        for sub_index in useful_sub_indexes:
            
            print sub_index,new_index
            
            sub_mask = (HO_sub_data == sub_index + 1)
            
            print np.sum(sub_mask == 1)
            
            
            full_HO_data[sub_mask] = new_index
            
            new_index = new_index + 1
            
        for cortl_index in useful_cortl_indexes:
            
            print cortl_index,new_index
            
            cortl_mask = (HO_cortl_data == cortl_index + 1)
            
            print np.sum(cortl_mask == 1)
            
            
            full_HO_data[cortl_mask] = new_index
            
            new_index = new_index + 1
            
        full_HO_data = np.array(full_HO_data,dtype = 'int')
        
        print "Original HO template shape:"
        print full_HO_data.shape, np.min(full_HO_data),np.max(full_HO_data)
        
        nib.save(nib.Nifti1Image(dataobj = full_HO_data,header = HO_sub_img_header,affine = HO_sub_img_affine),full_HO_img_file)
        
    else:
        
        full_HO_data = nib.load(full_HO_img_file).get_data()
        
    print np.unique(full_HO_data)
     
    #### reslicing to ref image
    resliced_full_HO_img_file = os.path.join(ROI_dir, "rHarvard-Oxford-cortl-sub-recombined-111regions.nii")
    
    print resliced_full_HO_img_file
    
    if not os.path.isfile(resliced_full_HO_img_file):
    
        reslice_full_HO = spm.Reslice()
        reslice_full_HO.inputs.in_file = full_HO_img_file
        reslice_full_HO.inputs.space_defining = img_file
        reslice_full_HO.inputs.out_file = resliced_full_HO_img_file
        #reslice_full_HO.inputs.interp = 0
        
        reslice_full_HO.run()
        
    ### loading results
    resliced_full_HO_img = nib.load(resliced_full_HO_img_file)
    
    resliced_full_HO_data = resliced_full_HO_img.get_data()
    print resliced_full_HO_data
    
    resliced_full_HO_data[np.isnan(resliced_full_HO_data)] = 0.0
    
    print np.unique(resliced_full_HO_data)
    
    
    
    resliced_full_HO_data = np.array(resliced_full_HO_data,dtype = int)
    
    print "Resliced HO template shape:"
    #print resliced_full_HO_data.shape
    
    #print resliced_full_HO_data
    
    print np.unique(resliced_full_HO_data)
    
    #### reading HO labels and concatenate
    
    ### sub
    HO_sub_labels_file = os.path.join(HO_dir,"HarvardOxford-Subcortical.xml")
    
    xmldoc_sub = minidom.parse(HO_sub_labels_file)
    
    HO_sub_labels = [s.firstChild.data for i,s in enumerate(xmldoc_sub.getElementsByTagName('label')) if i in useful_sub_indexes]
    
    print HO_sub_labels
    
    ### cortl
    HO_cortl_labels_file = os.path.join(HO_dir,"HarvardOxford-Cortical-Lateralized.xml")
    
    xmldoc_cortl = minidom.parse(HO_cortl_labels_file)
    
    HO_cortl_labels = [s.firstChild.data for s in xmldoc_cortl.getElementsByTagName('label') if i in useful_cortl_indexes]
    
    print HO_cortl_labels
    
    HO_labels = HO_sub_labels + HO_cortl_labels
    
    print len(HO_labels)
    
    HO_abbrev_labels = []
    
    for label in HO_labels:
    
        abbrev_label = ""
        
        split_label_parts = label.split(",")
        
        #print len(split_label_parts)
        
        for i_part,label_part in enumerate(split_label_parts):
        
            split_label = label_part.split()
            
            if i_part == 0:
                
                #print split_label
                
                ### left right
                if len(split_label) > 0:
                    
                    if split_label[0] == "Left":
                    
                        abbrev_label = "L."
                        
                    elif split_label[0] == "Right":
                    
                        abbrev_label = "R."
                        
                    else:
                        #print split_label[0]
                        
                        abbrev_label = split_label[0].title() 
                        #continue
                        
                        
                if len(split_label) > 1:
                
                    #print split_label[1]
                    
                    abbrev_label = abbrev_label + split_label[1][:5].title() 
                    
                if len(split_label) > 2:
                
                    #print split_label[2]
                    
                    abbrev_label = abbrev_label + split_label[2][:3].title() 
                    
                #if len(split_label) > 3 and split_label[3] != "":
                    
                    #for i in range(3,len(split_label)):
                    
                        #if split_label[i] != "":
                        
                            #print i,split_label[i]
                            
                        
            if i_part == 1:
                    
                split_label = split_label_parts[1].split()
                
                #print split_label
                
                ### left right
                if len(split_label) > 0 and split_label[0] != "":
                    
                    abbrev_label = abbrev_label + "." + split_label[0][:4].title() 
                    
                    #print abbrev_label
                
                if len(split_label) > 1:
                
                    #print split_label[1]
                    
                    abbrev_label = abbrev_label + split_label[1][:3].title() 
                    #0/0

        print abbrev_label
        
        HO_abbrev_labels.append(abbrev_label)
        #0/0
            
    print HO_abbrev_labels
    
    
    #0/0
    np_HO_abbrev_labels = np.array(HO_abbrev_labels,dtype = 'string')
    
    np_HO_labels = np.array(HO_labels,dtype = 'string')
    
    print np.unique(resliced_full_HO_data)
    
    template_indexes = np.unique(resliced_full_HO_data)[1:]
    
    print template_indexes
        
    print np_HO_labels.shape,np_HO_abbrev_labels.shape,template_indexes.shape
    
    #info_rois = np.hstack((np.unique(indexed_mask_rois_data)[1:].reshape(len(label_rois),1),np_full_label_rois,np_label_rois,rois_MNI_coords))
    #info_rois = np.hstack((np.unique(indexed_mask_rois_data)[1:].reshape(len(label_rois),1),rois_MNI_coords))
    info_template = np.hstack((template_indexes.reshape(len(HO_labels),1),np_HO_labels.reshape(len(HO_labels),1),np_HO_abbrev_labels.reshape(len(HO_labels),1)))
    #,rois_MNI_coords))
    
    print info_template
   
    info_template_file  =  os.path.join(ROI_dir, "info-Harvard-Oxford-reorg.txt")
    
    np.savetxt(info_template_file,info_template, fmt = '%s %s %s')
    
    return resliced_full_HO_img_file,info_template_file
    
        
def compute_labelled_mask_from_HO_sub(resliced_full_HO_img_file,info_template_file,export_dir):


    labeled_mask = nib.load(resliced_full_HO_img_file)
    
    print labeled_mask
    
    labeled_mask_data = labeled_mask.get_data()
    
    labeled_mask_header = labeled_mask.get_header().copy()
    
    labeled_mask_affine = np.copy(labeled_mask.get_affine())
    
    
    useful_indexes = range(2,6) + range(7,10) + range(11,18)
    
    print useful_indexes
    
    ROI_mask_data = np.zeros(shape = labeled_mask_data.shape,dtype = labeled_mask_data.dtype)
    
    for label_index in useful_indexes:
        
        roi_mask = (labeled_mask_data == label_index)
        
        print np.sum(roi_mask == True)
        
        ROI_mask_data[roi_mask] = label_index
                
    print np.unique(ROI_mask_data)
    
    
    if not os.path.exists(export_dir):
        os.makedirs(export_dir)
        
        
    ROI_mask_file = os.path.join(export_dir,"ROI_mask.nii")
    
    nib.save(nib.Nifti1Image(np.array(ROI_mask_data,dtype = int),labeled_mask_affine,labeled_mask_header),ROI_mask_file)

    
    labels = [line.strip().split(' ')[-1] for line in open(info_template_file)]
    
    np_labels = np.array(labels,dtype = 'string')
    
    print labels
    
    useful_labels = np_labels[np.array(useful_indexes,dtype = 'int')-1]
    
    print useful_labels
    
    ROI_mask_labels_file = os.path.join(export_dir,"ROI_mask_labels.txt")
    
    np.savetxt(ROI_mask_labels_file,useful_labels,fmt = "%s")
    
    
    ##export each ROI in a single file (mask)
    for label_index in useful_indexes:
    
        roi_mask = (labeled_mask_data == label_index)
    
        single_ROI_mask = nib.Nifti1Image(np.array(roi_mask,dtype = int),labeled_mask_affine,labeled_mask_header)
        
        print np_labels[label_index - 1]
        
        single_ROI_mask_file = os.path.join(export_dir,"ROI_mask_HO_" + np_labels[label_index - 1] + ".nii")
        
        nib.save(single_ROI_mask,single_ROI_mask_file)
        
    return ROI_mask_file,ROI_mask_labels_file
    
def segment_atlas_in_cubes(ROI_dir,ROI_cube_size,min_nb_voxels_in_neigh):
    
    from neuropype_graph.peak_labelled_mask import return_indexed_mask_neigh_within_binary_template
    
    resliced_atlas_file = os.path.join(ROI_dir,"rHarvard-Oxford-cortl-sub-recombined-111regions.nii")
    
    print resliced_atlas_file
    
    resliced_atlas_labels_file = os.path.join(ROI_dir,"info-Harvard-Oxford-reorg.txt")
    
    if not os.path.exists(resliced_atlas_file) or not os.path.exists(resliced_atlas_labels_file):
        
        print "Warning, missing atals: " + resliced_atlas_file
        
        return
        
    resliced_atlas = nib.load(resliced_atlas_file)
    
    resliced_atlas_data = resliced_atlas.get_data()
    
    reslice_atlas_header = resliced_atlas.get_header()
    
    reslice_atlas_affine = resliced_atlas.get_affine()
    
    
    print resliced_atlas_data.shape
    
    
    labels = [line.strip().split(' ')[-1] for line in open(resliced_atlas_labels_file)]
    
    np_labels = np.array(labels,dtype = 'string')
    
    print labels
    
    ############################ computing #####################################################
    
    list_selected_peaks_coords,indexed_mask_rois_data,label_rois = generate_continuous_mask_in_atlas(resliced_atlas_data,labels,ROI_cube_size,min_nb_voxels_in_neigh)

    print len(label_rois)
    
    template_indexes = np.array([resliced_atlas_data[coord[0],coord[1],coord[2]] for coord in list_selected_peaks_coords],dtype = 'int64')
    
    print template_indexes-1
    
    #label_rois = np_HO_abbrev_labels[template_indexes-1]
    full_label_rois = np_labels[template_indexes-1]
    
    #print label_rois2
    
    print label_rois
    
    #### exporting Rois image with different indexes 
    print np.unique(indexed_mask_rois_data)[1:].shape
    
    ROI_mask_prefix = "template_ROI_cube_size_" + str(ROI_cube_size)
    
    indexed_mask_rois_file = os.path.join(ROI_dir, "indexed_mask-" + ROI_mask_prefix + ".nii")
    
    nib.save(nib.Nifti1Image(dataobj = indexed_mask_rois_data,header = reslice_atlas_header,affine = reslice_atlas_affine),indexed_mask_rois_file)
    
    #### saving ROI coords as textfile
    
    coord_rois_file =  os.path.join(ROI_dir, "coords-" + ROI_mask_prefix + ".txt")

    np.savetxt(coord_rois_file,np.array(list_selected_peaks_coords,dtype = int), fmt = '%d')
    
    
    #### saving labels
    label_rois_file = os.path.join(ROI_dir,"labels-" + ROI_mask_prefix + ".txt")
    
    np_full_label_rois = np.array(full_label_rois,dtype = 'string').reshape(len(full_label_rois),1)
    
    print np_full_label_rois.shape
    
    np.savetxt(label_rois_file,np_full_label_rois, fmt = '%s')
    
    #info_rois = np.hstack((np.unique(indexed_mask_rois_data)[1:].reshape(len(label_rois),1),np_full_label_rois,np_label_rois,rois_MNI_coords))
    #info_rois = np.hstack((np.unique(indexed_mask_rois_data)[1:].reshape(len(label_rois),1),rois_MNI_coords))
    info_rois = np.hstack((np.unique(indexed_mask_rois_data)[1:].reshape(len(full_label_rois),1),np.array(list_selected_peaks_coords,dtype = int),np_full_label_rois))
    
    print info_rois
    info_rois_file = os.path.join(ROI_dir,"info-" + ROI_mask_prefix + ".txt")
   
    np.savetxt(info_rois_file,info_rois, fmt = '%s %s %s %s %s')
    
    
    return indexed_mask_rois_file,coord_rois_file
    
def generate_peaks(ROI_cube_size,mask_shape):
    
    from itertools import product
    
    ijk_list = [range(0,shape_i,ROI_cube_size) for shape_i in mask_shape]
    
    print ijk_list
    
    
    list_orig_peaks = [[i,j,k] for i,j,k in product(*ijk_list)]
    
    print list_orig_peaks
    
    return list_orig_peaks

    
def generate_continuous_mask_in_atlas(template_data,template_labels,ROI_cube_size,min_nb_voxels_in_neigh):
    
    from neuropype_graph.utils_dtype_coord import convert_np_coords_to_coords_dt
    from neuropype_graph.peak_labelled_mask import return_voxels_within_same_region

    template_data_shape = template_data.shape
    
    indexed_mask_rois_data = np.zeros(template_data_shape,dtype = 'int64') -1
    
    print indexed_mask_rois_data.shape
    
    list_orig_peak_coords = generate_peaks(ROI_cube_size,template_data_shape)
    
    label_rois = []
    
    list_selected_peaks_coords = []
    
    for orig_peak_coord in list_orig_peak_coords:
        
        print orig_peak_coord
        
        orig_peak_coord_np = np.array(orig_peak_coord)
        
        print orig_peak_coord_np
        
        list_voxel_coords,peak_template_roi_index = return_voxels_within_same_region(orig_peak_coord_np,ROI_cube_size,template_data,min_nb_voxels_in_neigh)
        
        print peak_template_roi_index
        
        if peak_template_roi_index > 0:
            
            neigh_coords = np.array(list_voxel_coords,dtype = 'int16')
    
            indexed_mask_rois_data[neigh_coords[:,0],neigh_coords[:,1],neigh_coords[:,2]] = len(list_selected_peaks_coords)
        
            label_rois.append(template_labels[peak_template_roi_index-1])
            
            list_selected_peaks_coords.append(orig_peak_coord_np)
                
            print len(list_selected_peaks_coords)
        
    return list_selected_peaks_coords,indexed_mask_rois_data,label_rois
    
