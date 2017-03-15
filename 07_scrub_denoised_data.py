__author__ = 'kanaan 18.12.2015'


import os
from utilities.utils import mkdir_path
from variables.subject_list import *
import subprocess
import shutil
import numpy as np


def scrub_data(population, workspace_dir):

    scub_subjects = []

    for subject in population:
        print '###############################################################################'
        print 'Scrubbing Denoised Data for subject %s' %subject
        print ''

        #input
        subject_dir       = os.path.join(workspace_dir, subject)

        native_residual_compor   = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_NATIVE_detrend_compcor_friston_bp_fwhm.nii.gz')
        native_residual_wmcsf    = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_NATIVE_detrend_wmcsf_friston_bp_fwhm.nii.gz')
        native_residual_global   = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_NATIVE_detrend_global_wmcsf_friston_bp_fwhm.nii.gz')

        mni_residual_compor   = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm.nii.gz')
        mni_residual_wmcsf    = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm.nii.gz')
        mni_residual_global   = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm.nii.gz')

        mni_aroma_compcor     = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp.nii.gz')
        mni_aroma_wmcsf       = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp.nii.gz')
        mni_aroma_global      = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp.nii.gz')

        #output
        scrub_dir =  os.path.join(subject_dir,'FUNC_DENOISE_SCRUB')
        mkdir_path(scrub_dir)
        os.chdir(scrub_dir)

        # scrubbing function
        def scrub(denoised_img, denoised_name):
            print '..Scrubbing Denoised image'
            print '----->', denoised_name[9:]
            scrubbed   = os.path.join(scrub_dir , '%s_scrubbed.nii.gz'%denoised_name)
            if not os.path.isfile(scrubbed):
                os.system("3dcalc -a %s%s -expr 'a' -prefix %s" %(denoised_img, in_frames_string, scrubbed))

        # get good frames amd limit to 150
        in_frames = []
        fd1d        =   np.genfromtxt(os.path.join(subject_dir,'QUALITY_CONTROL/FD.1D'))
        for frame, fd in enumerate(fd1d):
            if fd < 0.2:
                in_frames.append(frame)
        if len(in_frames) > 150:
            print '..Subject has more than 150 Good Frames (3.6 mins)'
            in_frames_string =  str(in_frames[0:150]).replace(" ","")

            # run scrubbing
            scrub(native_residual_compor, 'RESIDUAL_NATIVE_detrend_compcor_friston_bp_fwhm')
            scrub(native_residual_wmcsf,  'RESIDUAL_NATIVE_detrend_wmcsf_friston_bp_fwhm')
            scrub(native_residual_global, 'RESIDUAL_NATIVE_detrend_global_wmcsf_friston_bp_fwhm')

            scrub(mni_residual_compor, 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm')
            scrub(mni_residual_wmcsf,  'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm')
            scrub(mni_residual_global, 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm')

            scrub(mni_aroma_compcor, 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp')
            scrub(mni_aroma_wmcsf, 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp')
            scrub(mni_aroma_global, 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp')


        else:
            print 'Subject has less than 150 Good frames'

#scrub_data(['BM8X'], workspace_a)
#scrub_data(controls_a, workspace_a)
#scrub_data(patients_a, workspace_a)
scrub_data(patients_b, workspace_b)