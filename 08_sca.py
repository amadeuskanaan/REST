__author__ = 'kanaan 05.01.2016'

import os
from utilities.utils import *
from variables.subject_list import *

def run_seed_based_correlation(population, workspace_dir):

    print 'Running Seed Based Functional connectivity'

    count = 0
    for subject in population:
        count +=1

        print '###########################################'
        print '%s.Running Subject %s'%(count,subject)
        print ''

        subject_dir    = os.path.join(workspace_dir, subject)
        sca_dir        = os.path.join(subject_dir, 'SCA')
        mkdir_path(sca_dir)
        os.chdir(sca_dir)

        residual_compor   = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm.nii.gz')
        residual_wmcsf    = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm.nii.gz')
        residual_global   = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm.nii.gz')
        aroma_compcor     = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp.nii.gz')
        aroma_wmcsf       = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp.nii.gz')
        aroma_global      = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp.nii.gz')

        scrub_residual_compor   = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm_scrubbed.nii.gz')
        scrub_residual_wmcsf    = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm_scrubbed.nii.gz')
        scrub_residual_global   = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm_scrubbed.nii.gz')
        scrub_aroma_compcor     = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp_scrubbed.nii.gz')
        scrub_aroma_wmcsf       = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp_scrubbed.nii.gz')
        scrub_aroma_global      = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp_scrubbed.nii.gz')

        putamen  = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_PUTAMEN_L.nii.gz')
        caudate  = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_CAUDATE_L.nii.gz')
        pallidum = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_PALLIDUM_L.nii.gz')
        thalamus = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_THALAMUS_L.nii.gz')
        first    = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_FIRST.nii.gz')


        def run_sca(pproc, mask, mask_name, string):
            sca_dir        = os.path.join(subject_dir, 'SCA/%s'%mask_name)
            mkdir_path(sca_dir)
            os.chdir(sca_dir)

            print string
            if not os.path.isfile('SCA_Z_%s.nii.gz'%string):

                ##'1. Extract Timeseries'
                os.system('3dROIstats -quiet -mask_f2short -mask %s %s > rest_native.1d'%(mask, pproc))

                ##'2. Compute voxel-wise correlation with Seed Timeseries'
                os.system('3dfim+ -input %s -ideal_file rest_native.1d -out Correlation -bucket corr.nii.gz'%pproc)

                ##'3. Z-transform correlation'
                eq = ['log((a+1)/(a-1))/2']
                os.system('3dcalc -a corr.nii.gz -expr %s -prefix SCA_Z_%s.nii.gz'%(eq,string))
                os.system('rm -rf rest_native.1d corr.nii.gz')


        run_sca(residual_compor, putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm')
        run_sca(residual_wmcsf,  putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm')
        run_sca(residual_global, putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm')
        run_sca(aroma_compcor, putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp')
        run_sca(aroma_wmcsf,  putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp')
        run_sca(aroma_global, putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp')
        run_sca(scrub_residual_compor, putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm_scrubbed')
        run_sca(scrub_residual_wmcsf, putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm_scrubbed')
        run_sca(scrub_residual_global, putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm_scrubbed')
        run_sca(scrub_aroma_compcor, putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp_scrubbed')
        run_sca(scrub_aroma_global, putamen, 'PUTAMEN', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp_scrubbed')

        run_sca(residual_compor, caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm')
        run_sca(residual_wmcsf,  caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm')
        run_sca(residual_global, caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm')
        run_sca(aroma_compcor, caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp')
        run_sca(aroma_wmcsf,  caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp')
        run_sca(aroma_global, caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp')
        run_sca(scrub_residual_compor, caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm_scrubbed')
        run_sca(scrub_residual_wmcsf, caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm_scrubbed')
        run_sca(scrub_residual_global, caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm_scrubbed')
        run_sca(scrub_aroma_compcor, caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp_scrubbed')
        run_sca(scrub_aroma_global, caudate, 'CAUDATE', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp_scrubbed')

        run_sca(residual_compor, pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm')
        run_sca(residual_wmcsf,  pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm')
        run_sca(residual_global, pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm')
        run_sca(aroma_compcor, pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp')
        run_sca(aroma_wmcsf,  pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp')
        run_sca(aroma_global, pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp')
        run_sca(scrub_residual_compor, pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm_scrubbed')
        run_sca(scrub_residual_wmcsf, pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm_scrubbed')
        run_sca(scrub_residual_global, pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm_scrubbed')
        run_sca(scrub_aroma_compcor, pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp_scrubbed')
        run_sca(scrub_aroma_global, pallidum, 'PALLIDUM', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp_scrubbed')

        run_sca(residual_compor, thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm')
        run_sca(residual_wmcsf,  thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm')
        run_sca(residual_global, thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm')
        run_sca(aroma_compcor, thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp')
        run_sca(aroma_wmcsf,  thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp')
        run_sca(aroma_global, thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp')
        run_sca(scrub_residual_compor, thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm_scrubbed')
        run_sca(scrub_residual_wmcsf, thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm_scrubbed')
        run_sca(scrub_residual_global, thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm_scrubbed')
        run_sca(scrub_aroma_compcor, thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp_scrubbed')
        run_sca(scrub_aroma_global, thalamus, 'THALAMUS', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp_scrubbed')

        run_sca(residual_compor, first, 'FIRST', 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm')
        run_sca(residual_wmcsf,  first, 'FIRST', 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm')
        run_sca(residual_global, first, 'FIRST', 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm')
        run_sca(aroma_compcor, first, 'FIRST', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp')
        run_sca(aroma_wmcsf,  first, 'FIRST', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp')
        run_sca(aroma_global, first, 'FIRST', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp')
        run_sca(scrub_residual_compor, first, 'FIRST', 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm_scrubbed')
        run_sca(scrub_residual_wmcsf, first, 'FIRST', 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm_scrubbed')
        run_sca(scrub_residual_global, first, 'FIRST', 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm_scrubbed')
        run_sca(scrub_aroma_compcor, first, 'FIRST', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp_scrubbed')
        run_sca(scrub_aroma_global, first, 'FIRST', 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp_scrubbed')


# run_seed_based_correlation(controls_a, workspace_a)
# run_seed_based_correlation(patients_a, workspace_a)
run_seed_based_correlation(patients_b, workspace_b)