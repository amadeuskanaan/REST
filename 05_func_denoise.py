__author__ = 'kanaan 03.12.2015'

import os
from variables.subject_list import *
from utilities.utils import *
from utilities.nuisance import *
from utilities.bandpass import *
from utilities.ica_aroma import *
import shutil
import nipype.interfaces.ants as ants


def run_functional_denoising(population, datadir, workspace_dir):

    for subject in population:
        print '###############################################################################'
        print 'Denoising Functional Data for subject %s' %subject
        print ''

        #input
        subject_dir       = os.path.join(workspace_dir, subject)
        func_native       = os.path.join(subject_dir, 'FUNC_PPROC/REST_PPROC_NATIVE_BRAIN.nii.gz')
        func_native_mask  = os.path.join(subject_dir, 'FUNC_PPROC/REST_PPROC_NATIVE_BRAIN_mask_ero.nii.gz')
        func_native_mean  = os.path.join(subject_dir, 'FUNC_PPROC/REST_PPROC_NATIVE_BRAIN_mean.nii.gz')
        func_native_gm    = os.path.join(subject_dir, 'FUNC_TRANSFORM/NATIVE_FUNC_GM.nii.gz')
        func_native_wm    = os.path.join(subject_dir, 'FUNC_TRANSFORM/NATIVE_FUNC_WM.nii.gz')
        func_native_csf   = os.path.join(subject_dir, 'FUNC_TRANSFORM/NATIVE_FUNC_CSF.nii.gz')

        func_mni          = os.path.join(subject_dir, 'FUNC_TRANSFORM/REST_PPROC_MNI2mm_BRAIN.nii.gz')
        func_mni_mask     = os.path.join(subject_dir, 'FUNC_TRANSFORM/REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz')
        func_mni_mean     = os.path.join(subject_dir, 'FUNC_TRANSFORM/REST_PPROC_MNI2mm_BRAIN_mean.nii.gz')
        func_mni_gm       = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_GM.nii.gz')
        func_mni_wm       = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_WM.nii.gz')
        func_mni_csf      = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_CSF.nii.gz')

        anat_native   = os.path.join(subject_dir, 'ANATOMICAL/MP2RAGE.nii.gz')
        anat2func_aff = os.path.join(subject_dir, 'FUNC_TRANSFORM/FUNC2ANAT/ANAT2FUNC_XFM.mat')
        anat2mni_aff  = os.path.join(subject_dir, 'ANATOMICAL/transform0Affine.mat')
        anat2mni_wrp  = os.path.join(subject_dir, 'ANATOMICAL/transform1InverseWarp.nii.gz')



        #######################################################################################
        #output
        denoise_dir =  os.path.join(subject_dir, 'FUNC_DENOISE')
        mkdir_path(denoise_dir)
        os.chdir(denoise_dir)

        #  calculate friston parameters
        movpar  = os.path.join(subject_dir,'FUNC_PPROC/MOTION_CORRECTION/REST_DISCO_MOCO_2.1D')
        friston = os.path.join(subject_dir, 'FUNC_DENOISE/FRISTON_24.1D')
        #print '.... Calculating Friston-24 movpar'
        if not os.path.isfile(friston):
            calc_friston_twenty_four(movpar)

        #######################################################################################
        'NATIVE FUNC NUISANCE'
        #Extract Tissue Data
        print 'A. Native Func Nuisance Signal Regression'
        print '... Extracting Tissue mask Signals from native brain'
        print ''
        native_tissuesig_dir =  os.path.join(subject_dir, 'FUNC_DENOISE/TISSUE_SIGNALS_NATIVE')
        mkdir_path(native_tissuesig_dir)
        os.chdir(native_tissuesig_dir)

        #Erode WM
        os.system('fslmaths %s -ero native_WM_ero.nii.gz'%func_native_wm)
        func_native_wm = os.path.join(native_tissuesig_dir, 'native_WM_ero.nii.gz')

        native_wm_sig  = os.path.join(native_tissuesig_dir, 'NUISANCE_SIGNALS_WM.npy')
        native_gm_sig  = os.path.join(native_tissuesig_dir, 'NUISANCE_SIGNALS_GM.npy')
        native_csf_sig = os.path.join(native_tissuesig_dir, 'NUISANCE_SIGNALS_CSF.npy')


        # LIMIT CSF to ventricles
        if not os.path.isfile(native_gm_sig):
            native_lateral_venctivles = os.path.join(native_tissuesig_dir, 'HO_CSF_NATIVE_BIN.nii.gz')
            if not os.path.isfile(native_lateral_venctivles):
                #warp tp native space
                os.system('WarpImageMultiTransform 3 %s HO_CSF_ANAT.nii.gz -R %s -i %s %s'
                          %(mni_lateral_ventricles, anat_native, anat2mni_aff, anat2mni_wrp))
                print func_native_mean
                print  anat2func_aff
                os.system('flirt -in HO_CSF_ANAT.nii.gz -ref %s -applyxfm -init %s -out HO_CSF_NATIVE'
                          %(func_native_mean, anat2func_aff))
                os.system('fslmaths HO_CSF_NATIVE.nii.gz -thr 0.5 -bin  HO_CSF_NATIVE_BIN')

            extract_tissue_signals(data_file            = func_native,
                                   ventricles_mask_file = native_lateral_venctivles,
                                   wm_seg_file          = func_native_wm,
                                   csf_seg_file         = func_native_csf,
                                   gm_seg_file          = func_native_gm,
                                   wm_threshold         =0.01,
                                   csf_threshold        =0.01,
                                   gm_threshold         =0.0)

        ####################
        # 1a NATIVE  Compcor
        print '-A1. Calculating Native Func Residuals [ Compcor/Friston/Linear/Quadratic ]'

        NATIVE_COMPCOR_DIR =  os.path.join(denoise_dir, 'NUISANCE_NATIVE_COMPCOR')
        mkdir_path(NATIVE_COMPCOR_DIR)
        os.chdir(NATIVE_COMPCOR_DIR)

        if not os.path.isfile('./residual.nii.gz'):
            calc_residuals(subject  = func_native, selector  = nuisance_compcor,
                           wm_sig_file = native_wm_sig, csf_sig_file = native_csf_sig, gm_sig_file = native_gm_sig,
                           motion_file = friston, compcor_ncomponents = 5)
        print ' Filtering and Smoothing FWHM=4'
        if not os.path.isfile('../RESIDUAL_NATIVE_detrend_compcor_friston_bp_fwhm.nii.gz'):
            #bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.001,0.01), sample_period = None)
            bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.008,0.1), sample_period = None)

            # Smooth data
            os.system('fslmaths bandpassed_demeaned_filtered.nii.gz -kernel gauss 1.698644 -fmean ../RESIDUAL_NATIVE_detrend_compcor_friston_bp_fwhm.nii.gz')

        ##################
        # 1b NATIVE WMCSF
        print '-A2. Calculating Native Func Residuals [ WMCSF/Friston/Linear/Quadratic ]'
	print ''
        NATIVE_WMCSF_DIR =  os.path.join(denoise_dir, 'NUISANCE_NATIVE_WMCSF')
        mkdir_path(NATIVE_WMCSF_DIR)
        os.chdir(NATIVE_WMCSF_DIR)

        if not os.path.isfile('./residual.nii.gz'):
            calc_residuals(subject  = func_native, selector = nuisance_wmcsf,
                           wm_sig_file = native_wm_sig, csf_sig_file = native_csf_sig, gm_sig_file  = native_gm_sig,
                           motion_file  = friston, compcor_ncomponents = 5)
        print '     Filtering and Smoothing FWHM=4'
        if not os.path.isfile('../RESIDUAL_NATIVE_detrend_wmcsf_friston_bp_fwhm.nii.gz'):
            #bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.001,0.01), sample_period = None)
            bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.008,0.1), sample_period = None)
            # Smooth data
            os.system('fslmaths bandpassed_demeaned_filtered.nii.gz -kernel gauss 1.698644 -fmean ../RESIDUAL_NATIVE_detrend_wmcsf_friston_bp_fwhm.nii.gz')

        ##################
        # 1c NATIVE GLOBAL
        print '-A3. Calculating Native Func Residuals [ GLOBAL/WMCSF/Friston/Linear/Quadratic ]'
	print ''
        NATIVE_GLOBAL_DIR =  os.path.join(denoise_dir, 'NUISANCE_NATIVE_GLOBAL')
        mkdir_path(NATIVE_GLOBAL_DIR)
        os.chdir(NATIVE_GLOBAL_DIR)

        if not os.path.isfile('./residual.nii.gz'):
            calc_residuals(subject  = func_native, selector = nuisance_global,
                           wm_sig_file = native_wm_sig, csf_sig_file = native_csf_sig, gm_sig_file  = native_gm_sig,
                           motion_file  = friston, compcor_ncomponents = 5)
        print '     Filtering and Smoothing FWHM=4'
        if not os.path.isfile('../RESIDUAL_NATIVE_detrend_global_wmcsf_friston_bp_fwhm.nii.gz'):
            #bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.001,0.01), sample_period = None)
            bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.008,0.1), sample_period = None)
            # Smooth data
            os.system('fslmaths bandpassed_demeaned_filtered.nii.gz -kernel gauss 1.698644 -fmean ../RESIDUAL_NATIVE_detrend_global_wmcsf_friston_bp_fwhm.nii.gz')

        ################################################################################################################
        ################################################################################################################

        'MNI2mm FUNC NUISANCE'

        #  Extract Tissue Data
        print ''
        print 'B. MNI Func Nuisance Signal Regression'
        print '...Extracting Tissue mask Signals from MNI Func brain'
        TISSUESIGNAL_DIR_MNI =  os.path.join(subject_dir, 'FUNC_DENOISE/TISSUE_SIGNALS_MNI')
        mkdir_path(TISSUESIGNAL_DIR_MNI)
        os.chdir(TISSUESIGNAL_DIR_MNI)

        #Erode WM
        os.system('fslmaths %s -ero mni_WM_ero.nii.gz'%func_mni_wm)
        func_native_wm = os.path.join(TISSUESIGNAL_DIR_MNI, 'mni_WM_ero.nii.gz')

        mni_wm_sig  = os.path.join(TISSUESIGNAL_DIR_MNI, 'NUISANCE_SIGNALS_WM.npy')
        mni_gm_sig  = os.path.join(TISSUESIGNAL_DIR_MNI, 'NUISANCE_SIGNALS_GM.npy')
        mni_csf_sig = os.path.join(TISSUESIGNAL_DIR_MNI, 'NUISANCE_SIGNALS_CSF.npy')

        if not os.path.isfile(mni_gm_sig):
            extract_tissue_signals(data_file            = func_mni,
                                   ventricles_mask_file = mni_lateral_ventricles,
                                   wm_seg_file          = func_mni_wm,
                                   csf_seg_file         = func_mni_csf,
                                   gm_seg_file          = func_mni_gm,
                                   wm_threshold         =0.01,
                                   csf_threshold        =0.01,
                                   gm_threshold         =0.0)

        ####################
        #  MNI  Compcor
        print '-B1. Calculating MNI Func Residuals [ Compcor/Friston/Linear/Quadratic ]'
        print ''
        MNI_COMPCOR_DIR =  os.path.join(denoise_dir, 'NUISANCE_MNI_COMPCOR')
        mkdir_path(MNI_COMPCOR_DIR)
        os.chdir(MNI_COMPCOR_DIR)
        if not os.path.isfile('./residual.nii.gz'):
            calc_residuals(subject        = func_mni, selector       = nuisance_compcor,
                           wm_sig_file = mni_wm_sig, csf_sig_file= mni_csf_sig, gm_sig_file = mni_gm_sig,
                           motion_file = friston, compcor_ncomponents = 5)
        print '     Filtering and Smoothing FWHM=4'
        if not os.path.isfile('../RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm.nii.gz'):
            #bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.001,0.01), sample_period = None)
            bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.008,0.1), sample_period = None)
            # Smooth data
            os.system('fslmaths bandpassed_demeaned_filtered.nii.gz -kernel gauss 1.698644 -fmean ../RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm.nii.gz')

        ####################
        # MNI  WMCSF
        print '-B2. Calculating MNI Func Residuals [ WMCSF/Friston/Linear/Quadratic ]'
        print ''
        MNI_WMCSF_DIR =  os.path.join(denoise_dir, 'NUISANCE_MNI_WMCSF')
        mkdir_path(MNI_WMCSF_DIR)
        os.chdir(MNI_WMCSF_DIR)
        if not os.path.isfile('./residual.nii.gz'):
            calc_residuals(subject = func_mni, selector= nuisance_wmcsf,
                           wm_sig_file = mni_wm_sig, csf_sig_file= mni_csf_sig, gm_sig_file = mni_gm_sig,
                           motion_file = friston, compcor_ncomponents = 5)
        print '     Filtering and Smoothing FWHM=4'
        if not os.path.isfile('../RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm.nii.gz'):
            # bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.001,0.01), sample_period = None)
            bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.008,0.1), sample_period = None)
            # Smooth data
            os.system('fslmaths bandpassed_demeaned_filtered.nii.gz -kernel gauss 1.698644 -fmean ../RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm.nii.gz')

        ####################
        #  MNI  GLOBAL
        print '-B3. Calculating MNI Func Residuals [GLOBAL/WMCSF/Friston/Linear/Quadratic ]'
        print ''
        MNI_GLOBAL_DIR =  os.path.join(denoise_dir, 'NUISANCE_MNI_GLOBAL')
        mkdir_path(MNI_GLOBAL_DIR)
        os.chdir(MNI_GLOBAL_DIR)

        #Nuisance
        if not os.path.isfile('./residual.nii.gz'):
            calc_residuals(subject = func_mni, selector= nuisance_global,
                           wm_sig_file = mni_wm_sig, csf_sig_file= mni_csf_sig, gm_sig_file = mni_gm_sig,
                           motion_file = friston, compcor_ncomponents = 5)
        # BP filter
        print '     Filtering and Smoothing FWHM=4'
        if not os.path.isfile('../RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm.nii.gz'):
            # bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.001,0.01), sample_period = None)
            bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.008,0.1), sample_period = None)
            # Smooth data
            os.system('fslmaths bandpassed_demeaned_filtered.nii.gz -kernel gauss 1.698644 -fmean ../RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm.nii.gz')

        ################################################################################################################
        ################################################################################################################
        ################################################################################################################
        ################################################################################################################
        # MNI  FUNC AROMA NUISANCE

        print ''
        print 'C. Denoising with AROMA'
        print ''
        MNI_AROMA_CC_DIR =  os.path.join(denoise_dir, 'NUISANCE_MNI_AROMA_COMPCOR')
        mkdir_path(MNI_AROMA_CC_DIR)
        os.chdir(MNI_AROMA_CC_DIR)


        # run aroma on smoothed mni func#
        func_mni_aroma = os.path.join(MNI_AROMA_CC_DIR, 'denoised_func_data_nonaggr.nii.gz')
        print '     Running ICA'
        if not os.path.isfile(func_mni_aroma):
            os.system('fslmaths %s -kernel gauss 1.698644 -fmean REST_PPROC_MNI2mm_BRAIN_fwhm.nii.gz' %func_mni)

            ica_aroma_denoise(fslDir  =  '/usr/share/fsl/5.0/bin/',
                              inFile  =  './REST_PPROC_MNI2mm_BRAIN_fwhm.nii.gz',
                              mask    =  func_mni_mask,
                              dim     =  0,
                              TR      =  1.3999,
                              mc      =  movpar,
                              denType = 'nonaggr',
                              affmat  =  [],
                              warp    =  [])

        # ####################
         # MNI AROMA COMPCOR
        print '-C1. Calculating Residuals AROMA + [Compcor/Friston/Linear/Quadratic ]'
        print ''
        if not os.path.isfile('./residual.nii.gz'):
            calc_residuals(subject = func_mni_aroma, selector = nuisance_compcor,
                           wm_sig_file = mni_wm_sig, csf_sig_file= mni_csf_sig, gm_sig_file = mni_gm_sig,
                           motion_file = friston, compcor_ncomponents = 5)
        
        print '     Filtering and Smoothing FWHM=4'
        if not os.path.isfile('../RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp.nii.gz'):
            # bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.001,0.01), sample_period = None)
            bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.008,0.1), sample_period = None)
            # Smooth data
            os.system('cp bandpassed_demeaned_filtered.nii.gz ../RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp.nii.gz')
        

        ####################
        # MNI AROMA WMCSF
        MNI_AROMA_WMCSF_DIR =  os.path.join(denoise_dir, 'NUISANCE_MNI_AROMA_WMCSF')
        mkdir_path(MNI_AROMA_WMCSF_DIR)
        os.chdir(MNI_AROMA_WMCSF_DIR)

        print '-C2. Calculating Residuals AROMA + [WMCSF/Friston/Linear/Quadratic ]'
        print ''
        if not os.path.isfile('./residual.nii.gz'):
            calc_residuals(subject = func_mni_aroma, selector = nuisance_wmcsf,
                           wm_sig_file = mni_wm_sig, csf_sig_file= mni_csf_sig, gm_sig_file = mni_gm_sig,
                           motion_file = friston, compcor_ncomponents = 5)
        print '     Filtering and Smoothing FWHM=4'
        if not os.path.isfile('../RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp.nii.gz'):
            # bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.001,0.01), sample_period = None)
            bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.008,0.1), sample_period = None)
            # Smooth data
            os.system('cp bandpassed_demeaned_filtered.nii.gz ../RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp.nii.gz')

        ####################
        # MNI AROMA GLOBAL
        MNI_AROMA_GLOBAL_DIR =  os.path.join(denoise_dir, 'NUISANCE_MNI_AROMA_GLOBAL')
        mkdir_path(MNI_AROMA_GLOBAL_DIR)
        os.chdir(MNI_AROMA_GLOBAL_DIR)

        print '-C3. Calculating Residuals AROMA + [GLOBAL/WMCSF/Friston/Linear/Quadratic ]'
        print ''
        if not os.path.isfile('./residual.nii.gz'):
            calc_residuals(subject = func_mni_aroma, selector = nuisance_global,
                           wm_sig_file = mni_wm_sig, csf_sig_file= mni_csf_sig, gm_sig_file = mni_gm_sig,
                           motion_file = friston, compcor_ncomponents = 5)
        print '     Filtering and Smoothing FWHM=4'
        if not os.path.isfile('../RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp.nii.gz'):
            # bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.001,0.01), sample_period = None)
            bandpass_voxels(realigned_file= './residual.nii.gz', bandpass_freqs= (0.008,0.1), sample_period = None)
            # Smooth data
            os.system('cp bandpassed_demeaned_filtered.nii.gz ../RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp.nii.gz')

# run_functional_denoising(['HR8T'], controls_datadir_a, workspace_a)
# run_functional_denoising(controls_a, controls_datadir_a, workspace_a)
run_functional_denoising(patients_a, patients_datadir_a, workspace_a)
# run_functional_denoising(patients_b, patients_datadir_b, workspace_b)

