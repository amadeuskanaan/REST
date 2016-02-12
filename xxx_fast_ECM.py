__author__ = 'kanaan 11.12.2015'

import os
from variables.subject_list import *
from utilities.utils import mkdir_path
from utilities.nuisance import *
from utilities.bandpass import *
import shutil
import subprocess
import commands

def run_eigenvector_centrality(population, workspace_dir):

    for subject in population:
        print '###############################################################################'
        print 'Running Eigenvector Centrality Mapping subject %s' %subject
        print ''

        #input
        subject_dir       = os.path.join(workspace_dir, subject)
        residual_compor   = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm.nii.gz')
        residual_wmcsf    = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm.nii.gz')
        residual_global   = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm.nii.gz')

        aroma_compcor     = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp.nii.gz')
        aroma_wmcsf       = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp.nii.gz')
        aroma_global      = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp.nii.gz')

        scrubbed_residual_compor   = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm_scrubbed.nii.gz')
        scrubbed_residual_wmcsf    = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm_scrubbed.nii.gz')
        scrubbed_residual_global   = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm_scrubbed.nii.gz')

        scrubbed_aroma_compcor     = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp_scrubbed.nii.gz')
        scrubbed_aroma_wmcsf       = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp_scrubbed.nii.gz')
        scrubbed_aroma_global      = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp_scrubbed.nii.gz')


        group_gm_mask     = '/SCR4/workspace/project_GluRest/OUT_DIR_A/GluConnectivity/GM2mm_bin.nii'


        def run_fast_ecm(nuisance_file, gm_mask, string):
            if not os.path.isfile(os.path.join(subject_dir, 'FAST_ECM',string, 'zscore_fastECM.nii.gz')):
                #create output dir
                dir =  os.path.join(subject_dir,'FAST_ECM', string)
                mkdir_path(dir)
                os.chdir(dir)

                #copy nuisance file locally
                shutil.copy(nuisance_file, './residual.nii.gz')

                #gunzip
                if not os.path.isfile('residual.nii'):
                    os.system('gunzip residual.nii.gz')
                    os.system('rm -rf residual.nii.gz')

                # run Fast ECM
                pproc = os.path.join(dir, 'residual.nii')
                matlab_cmd = ['matlab',  '-version', '8.2', '-nodesktop' ,'-nosplash'  ,'-nojvm' ,'-r "fastECM(\'%s\', \'1\', \'1\', \'1\', \'20\', \'%s\') ; quit;"' %(pproc, gm_mask)]
                print '    ... Running ECM'
                subprocess.call(matlab_cmd)

                def z_score_centrality(image, outname):

                    print '    ... z-scoring %s'%outname
                    # zscore fastECM image
                    std  = commands.getoutput('fslstats %s -k %s -s | awk \'{print $1}\''%(image, group_gm_mask))
                    mean = commands.getoutput('fslstats %s -k %s -m | awk \'{print $1}\''%(image, group_gm_mask))
                    os.system('fslmaths %s -sub %s -div %s -mas %s %s'%(image, mean, std, group_gm_mask, outname))

                z_score_centrality('residual_fastECM.nii', 'zscore_fastECM')
                z_score_centrality('residual_degCM.nii'  , 'zscore_degCM')
                z_score_centrality('residual_normECM.nii', 'zscore_normECM')
                z_score_centrality('residual_rankECM.nii', 'zscore_rankECM')

        print '1. Nuisance COMPCOR ECM'
        run_fast_ecm(residual_compor, group_gm_mask, 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm')
        print '.'
        run_fast_ecm(scrubbed_residual_compor, group_gm_mask, 'RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm_scrubbed')
        print '.'

        #######################
        print '2. Nuisance WMCSF ECM'
        run_fast_ecm(residual_wmcsf, group_gm_mask, 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm')
        print '.'
        run_fast_ecm(scrubbed_residual_wmcsf, group_gm_mask, 'RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm_scrubbed')
        print '.'

        #######################
        print '3. Nuisance GLOBAL ECM'
        run_fast_ecm(residual_global, group_gm_mask, 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm')
        print '.'
        run_fast_ecm(scrubbed_residual_global, group_gm_mask, 'RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm_scrubbed')
        print '.'

        #######################
        print '4. Nuisance AROMA-COMPCOR ECM'
        run_fast_ecm(aroma_compcor, group_gm_mask, 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp')
        print '.'
        run_fast_ecm(scrubbed_aroma_compcor, group_gm_mask, 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp_scrubbed')
        print '.'

        #######################
        print '5. Nuisance AROMA-WMCSF ECM'
        run_fast_ecm(aroma_wmcsf, group_gm_mask, 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp')
        print '.'
        run_fast_ecm(scrubbed_aroma_wmcsf, group_gm_mask, 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp_scrubbed')
        print '.'

        #######################
        print '6. Nuisance AROMA-GLOBAL ECM'
        run_fast_ecm(aroma_global, group_gm_mask, 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp')
        print '.'
        run_fast_ecm(scrubbed_aroma_global, group_gm_mask, 'RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp_scrubbed')
        print '.'


# run_eigenvector_centrality(['BM8X'], workspace_a)
run_eigenvector_centrality(controls_a, workspace_a)
run_eigenvector_centrality(patients_a, workspace_a)
run_eigenvector_centrality(patients_b, workspace_b)

