__author__ = 'kanaan 03.12.2015'

import os
from variables.subject_list import *
from utilities.utils import *
import shutil
import nipype.interfaces.ants as ants
import sys

# assert len(sys.argv)== 2
# subject_index=int(sys.argv[1])

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
                                       Anatomical Pre-processing
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def run_anatomical_preprocessing(population, datadir, workspace_dir):

    for subject in population:
        #subject = population[subject_index]
        print '###############################################################################'
        print 'Processing Anatomical Data for subject %s' %subject
        print ''

        # input
        anat   =  os.path.join(datadir, subject, 'NIFTI', 'MP2RAGE_DESKULL_RPI.nii.gz')
        gm     =  os.path.join(datadir, subject, 'NIFTI', 'TISSUE_CLASS_1_GM_OPTIMIZED.nii.gz')
        wm     =  os.path.join(datadir, subject, 'NIFTI', 'TISSUE_CLASS_2_WM_OPTIMIZED.nii.gz')
        csf    =  os.path.join(datadir, subject, 'NIFTI', 'TISSUE_CLASS_3_CSF_OPTIMIZED.nii.gz')
        first4 =  os.path.join(datadir, subject, 'NIFTI', 'FIRST_4d.nii.gz')
        first3 =  os.path.join(datadir, subject, 'NIFTI', 'FIRST_3d.nii.gz')

        #output
        anat_dir =  os.path.join(workspace_dir, subject, 'ANATOMICAL')
        mkdir_path(anat_dir)
        os.chdir(anat_dir)

        # Grabbing and Reorienting anatomical image and masks
        if not os.path.isfile(os.path.join(anat_dir, 'SUBCORTICAL_3D.nii.gz')):
            #Erode anat
            os.system('fslmaths %s -ero MP2RAGE.nii.gz' %anat)

            # reorient masks
            os.system('fslswapdim %s RL PA IS GM.nii.gz' %gm)
            os.system('fslswapdim %s RL PA IS WM.nii.gz' %wm)
            os.system('fslswapdim %s RL PA IS CSF.nii.gz' %csf)
            os.system('fslswapdim %s RL PA IS SUBCORTICAL_4D.nii.gz' %first4)
            os.system('fslswapdim %s RL PA IS SUBCORTICAL_3D.nii.gz' %first3)


        # anat2mni
        print 'MP2RAGE to MNI Linear-Registration '
        if not os.path.isfile('./MP2RAGE_MNI1mm.nii.gz'):
            # os.system('flirt -in MP2RAGE.nii.gz -ref %s -out MP2RAGE_MNI1mm_Linear.nii.gz -omat MP2RAGE_MNI1mm_Linear.mat -cost mutualinfo -dof 12' %mni_brain_1mm)


            ants_anat2mni = ants.Registration()
            ants_anat2mni.inputs.moving_image               = os.path.join(anat_dir, 'MP2RAGE.nii.gz')
            ants_anat2mni.inputs.fixed_image                =   mni_brain_1mm
            ants_anat2mni.inputs.dimension                  = 3                                          # Dimesion of input (default is 3)
            ants_anat2mni.inputs.use_histogram_matching     = True                          # Match hists of images before reg
            ants_anat2mni.inputs.winsorize_lower_quantile   = 0.01                        # Winsorize data based on quantilies (lower  value)
            ants_anat2mni.inputs.winsorize_upper_quantile   = 0.99                        # Winsorize data based on quantilies (higher value)
            ants_anat2mni.inputs.metric                     = ['MI','CC']                                   # Image metric(s) to be used at each stage
            ants_anat2mni.inputs.metric_weight              = [1,1]                                  # Modulates the per-stage weighting of the corresponding metric
            ants_anat2mni.inputs.radius_or_number_of_bins   = [32,4]                      # Number of bins in each stage for the MI and Mattes metric
            ants_anat2mni.inputs.sampling_strategy          = ['Regular',None]                   # Sampling strategy to use for the metrics {None, Regular, or Random}
            ants_anat2mni.inputs.sampling_percentage        = [0.25,0.25,None]                 # Defines the sampling strategy
            ants_anat2mni.inputs.number_of_iterations       = [[300,200,100], [50,30,20]]   #[[1000,500,250,100],[1000,500,250,100], [100,100,70,20]]  # Determines the convergence
            ants_anat2mni.inputs.convergence_threshold      = [1e-8,1e-9]                    # Threshold compared to the slope of the line fitted in convergence
            ants_anat2mni.inputs.convergence_window_size    = [10,15]                      # Window size of convergence calculations
            ants_anat2mni.inputs.transforms                 = ['Affine','SyN']                          # Selection of registration options. See antsRegistration documentation
            ants_anat2mni.inputs.transform_parameters       = [[0.1],[0.1,3,0]]               # Selection of registration options. See antsRegistration documentation
            ants_anat2mni.inputs.transform_parameters       = [[0.1],[0.1,3,0]]               # Fine-tuning for the different registration options
            ants_anat2mni.inputs.shrink_factors             = [[4,2,1],[4,2,1]]                     #Specify the shrink factor for the virtual domain (typically the fixed image) at each level
            ants_anat2mni.inputs.smoothing_sigmas           = [[2,1,0],[2,1,0]]                   # Specify the sigma of gaussian smoothing at each level
            ants_anat2mni.inputs.output_warped_image        =  os.path.join(anat_dir, 'MP2RAGE_MNI1mm.nii.gz')
            ants_anat2mni.run()

            # Apply warp field
            apply_ants_warp  = ants.WarpImageMultiTransform()
            apply_ants_warp.inputs.input_image = os.path.join(anat_dir, 'MP2RAGE.nii.gz')
            apply_ants_warp.inputs.reference_image    = mni_brain_1mm
            apply_ants_warp.inputs.transformation_series = ['%s/transform1Warp.nii.gz'%anat_dir,'%s/transform0Affine.mat'%anat_dir]
            apply_ants_warp.run()
            os.system('mv MP2RAGE_wimt.nii.gz  MP2RAGE_MNI1mm.nii.gz')

        # Optimize Masks
        if not os.path.isfile('./GM_optimized.nii.gz'):
            # WARP ATAG ATLAS TO NATIVE SPACE
            os.system('WarpImageMultiTransform 3 %s ATAG_NATIVE.nii.gz -R %s -i transform0Affine.mat transform1InverseWarp.nii.gz'
                       %(atag_atlas, os.path.join(anat_dir, 'MP2RAGE.nii.gz') ))
            os.system('fslmaths ATAG_NATIVE.nii.gz -thr 0.1 -bin  ATAG_NATIVE_BIN')

            # subtract atag mask from WM and CSF and add to GM
            os.system('fslmaths WM.nii.gz -sub ATAG_NATIVE_BIN.nii.gz -bin WM_optimized.nii.gz')
            os.system('fslmaths CSF.nii.gz -sub ATAG_NATIVE_BIN.nii.gz -bin CSF_optimized.nii.gz')
            os.system('fslmaths GM.nii.gz -add ATAG_NATIVE_BIN.nii.gz -bin GM_optimized.nii.gz')

        # Separate FIRST Subcortical Masks

run_anatomical_preprocessing(['HCTT'], controls_datadir_a, workspace_a)
#run_anatomical_preprocessing(controls_a, controls_datadir_a, workspace_a)
#run_anatomical_preprocessing(patients_a, patients_datadir_a, workspace_a)
#run_anatomical_preprocessing(patients_b, patients_datadir_b, workspace_b)
