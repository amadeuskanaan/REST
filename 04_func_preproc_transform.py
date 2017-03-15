__author__ = 'kanaan 03.12.2015'

import os
from variables.subject_list import *
from utilities.utils import *
import shutil
import nipype.interfaces.ants as ants
from utilities.bandpass import  *



def run_functional_transform(population, datadir, workspace_dir, freesurfer_dir):

    for subject in population:
        print '###############################################################################'
        print 'Transforming Functional Data for subject %s' %subject
        print ''

        #input
        subject_dir     =  os.path.join(workspace_dir, subject)
        pproc_dir = os.path.join(subject_dir, 'FUNC_PPROC')
        func_pproc_mean = os.path.join(subject_dir, 'FUNC_PPROC/REST_PPROC_NATIVE_BRAIN_mean.nii.gz')
        func_pproc_mask = os.path.join(subject_dir, 'FUNC_PPROC/REST_PPROC_NATIVE_BRAIN_mask_ero.nii.gz')

        anat      = os.path.join(subject_dir, 'ANATOMICAL/MP2RAGE.nii.gz')
        gm        = os.path.join(subject_dir, 'ANATOMICAL/GM_optimized.nii.gz')
        wm        = os.path.join(subject_dir, 'ANATOMICAL/WM_optimized.nii.gz')
        csf       = os.path.join(subject_dir, 'ANATOMICAL/CSF_optimized.nii.gz')
        atag      = os.path.join(subject_dir, 'ANATOMICAL/ATAG_NATIVE_BIN.nii.gz')
        first_4d  = os.path.join(subject_dir, 'ANATOMICAL/SUBCORTICAL_4D.nii.gz')

        warp     =  os.path.join(subject_dir, 'ANATOMICAL/transform1Warp.nii.gz')
        affine   =  os.path.join(subject_dir, 'ANATOMICAL/transform0Affine.mat')
        warp_inv =  os.path.join(subject_dir, 'ANATOMICAL/transform1InverseWarp.nii.gz')

        func2nat_dir =  os.path.join(subject_dir, 'FUNC_TRANSFORM/FUNC2ANAT')
        transform_dir =  os.path.join(subject_dir, 'FUNC_TRANSFORM')
        mkdir_path(func2nat_dir)
        os.chdir(func2nat_dir)

         ############################################################################################################
        # FUNC_2.3mm --> ANAT_2mm

        # Calculate func -- > anat linear_xfm
        # runs flirt using a two step procedure (Mutual Information and Boundary based registration )

        print '....Transforming functional image to anatomical space'
        if not os.path.isfile('FUNC2ANAT_UNI.nii.gz'):
            # 1- resample anat to 2.3mm
            os.system('flirt -in %s -ref %s -out ANAT_2.3mm -applyisoxfm 2.3 -datatype float'%(anat,anat))

            # 2- register func to anat - step 1 mutualinfo
            print '........Step 1: Mutual Info'
            os.system('flirt -in %s -ref ANAT_2.3mm.nii.gz -cost mutualinfo -noresample -dof 6 -out FUNC2ANAT_STEP1_MI.nii.gz -omat FUNC2ANAT_STEP1_MI.mat'%func_pproc_mean)

            # 3- register func to anat -step 2 bbr with mutualinfo init
            print '........Step 2: BBR'
            os.system('flirt -in %s -ref ANAT_2.3mm.nii.gz -cost bbr -dof 6 -init FUNC2ANAT_STEP1_MI.mat -noresample -schedule /usr/share/fsl/5.0/etc/flirtsch/bbr.sch '
                      '-wmseg %s -out FUNC2ANAT_STEP_2_BBR -omat FUNC2ANAT_STEP_2_BBR.mat' %(func_pproc_mean, wm ))

            # 4- Unify moco and linear regisrtation affines
            mkdir_path(os.path.join(func2nat_dir, 'MATS_UNIFIED'))
            mats_dir = os.path.join(func2nat_dir, 'MATS_UNIFIED')
            os.chdir(mats_dir)

            #split diso volume
            os.system('fslsplit %s/REST_DISCO.nii.gz -t'%pproc_dir)

            #combine motion and anat affines
            print '........Concatenating linear and motion affines'
            for i in xrange(0,417):
                frame = '{0:0>4}'.format(i)
                os.system('convert_xfm -omat UNIMAT_%s.mat -concat ../FUNC2ANAT_STEP_2_BBR.mat %s/MOTION_CORRECTION/MATS/MAT_%s'%(frame, pproc_dir,frame))
                os.system('flirt -in vol%s.nii.gz -ref ../ANAT_2.3mm.nii.gz -applyxfm -init UNIMAT_%s.mat -out flirt_uni_%s.nii.gz'%(frame, frame,frame))

            os.system('fslmerge -t ../FUNC2ANAT_UNI.nii.gz flirt_uni_*')
            os.system('rm -rf vol* flirt_uni_*')

        ############################################################################################################
        # FUNC_2.3mm --> MNI_2mm
        func2mni_dir =  os.path.join(subject_dir, 'FUNC_TRANSFORM/FUNC2MNI')
        mkdir_path(func2mni_dir)
        os.chdir(func2mni_dir)

        print '....Transforming functional image to MNI space'
        # FUNC_2 --> MNI_2mm
        if not os.path.isfile('REST_PPROC_MNI2mm.nii'):

            os.system('WarpTimeSeriesImageMultiTransform 4 %s resting_wtsimt.nii -R %s %s %s'
                      %(os.path.join(subject_dir, 'FUNC_TRANSFORM/FUNC2ANAT/FUNC2ANAT_UNI.nii.gz'),
                        mni_brain_2mm, warp, affine))
            os.system('cp resting_wtsimt.nii REST_PPROC_MNI2mm.nii')

        # Intensity normalization and deskulling
        print '.....Normalizing intensity to Mode 1000 and deskulling'
        os.chdir(transform_dir)
        if not os.path.isfile('REST_PPROC_MNI2mm_BRAIN.nii.gz'):
            os.system('fslmaths FUNC2MNI/REST_PPROC_MNI2mm.nii -ing 1000 REST_PPROC_MNI2mm.nii.gz')
            os.system('bet REST_PPROC_MNI2mm.nii.gz REST_PPROC_MNI2mm_BRAIN.nii.gz -f 0.50 -F -m -t -g 0.00')
            os.system('3dTstat -mean -prefix REST_PPROC_MNI2mm_BRAIN_mean.nii.gz  REST_PPROC_MNI2mm_BRAIN.nii.gz ')
            os.system('fslmaths REST_PPROC_MNI2mm_BRAIN_mask.nii.gz -ero -ero REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz')

        ############### CREATING PPRC IMAGE READY FOR RESOURCE ALLOCATION INDEX COMPUTATION
        print '.....Creating pproc/warped image ready for ICA '
        # ica ready file
        if not os.path.isfile('ICAready_REST_PPROC_MNI2mm_fwhm_bp.nii.gz'):
            # Smoothing and bp
            os.system('fslmaths FUNC2MNI/REST_PPROC_MNI2mm.nii -kernel gauss 1.698644 -fmean ICAready_REST_PPROC_MNI2mm_fwhm.nii.gz')
            bandpass_voxels(realigned_file= 'ICAready_REST_PPROC_MNI2mm_fwhm.nii.gz', bandpass_freqs= (0.008,0.1), sample_period = None)
            os.system('cp bandpassed_demeaned_filtered.nii.gz ICAready_REST_PPROC_MNI2mm_fwhm_bp.nii.gz')


        ############################################################################################################
        # ANAT MASKS to FUNC NATIVE
        os.chdir(func2nat_dir)
        print '.....Transforming anatomical masks to functional space'
        #invert xfm
        if not os.path.isfile(os.path.join(subject_dir, 'FUNC_TRANSFORM/NATIVE_FUNC_ATAG.nii.gz')):
            os.system('convert_xfm -omat ANAT2FUNC_XFM.mat -inverse FUNC2ANAT_STEP_2_BBR.mat')
            os.system('flirt -in %s -ref %s -out GM_FUNC.nii.gz -applyxfm -init ANAT2FUNC_XFM.mat'%(gm, func_pproc_mean,))
            os.system('flirt -in %s -ref %s -out WM_FUNC.nii.gz -applyxfm -init ANAT2FUNC_XFM.mat'%(wm, func_pproc_mean,))
            os.system('flirt -in %s -ref %s -out CSF_FUNC.nii.gz -applyxfm -init ANAT2FUNC_XFM.mat'%(csf, func_pproc_mean,))
            os.system('flirt -in %s -ref %s -out ATAG_FUNC.nii.gz -applyxfm -init ANAT2FUNC_XFM.mat'%(atag, func_pproc_mean,))


            os.system('fslmaths GM_FUNC.nii.gz -thr 0.3 -mul %s -bin ../NATIVE_FUNC_GM.nii.gz'%func_pproc_mask)
            os.system('fslmaths WM_FUNC.nii.gz -thr 0.7 -mul %s -bin WM_FUNCxx.nii.gz'%func_pproc_mask)
            os.system('fslmaths CSF_FUNC.nii.gz -thr 0.7 -mul %s -bin ../NATIVE_FUNC_CSF.nii.gz'%func_pproc_mask)
            os.system('fslmaths ATAG_FUNC.nii.gz -thr 0.3 -bin ../NATIVE_FUNC_ATAG.nii.gz')
            os.system('fslmaths WM_FUNCxx.nii.gz -sub ../NATIVE_FUNC_ATAG.nii.gz -bin ../NATIVE_FUNC_WM.nii.gz')

        if not os.path.isfile(os.path.join(subject_dir, 'FUNC_TRANSFORM/NATIVE_FUNC_FIRST.nii.gz')):
            os.system('fslroi %s NACC_l.nii.gz 0 1'%first_4d)
            os.system('fslroi %s AMYGDALA_l.nii.gz 1 1'%first_4d)
            os.system('fslroi %s CAUDATE_l.nii.gz 2 1'%first_4d)
            os.system('fslroi %s HIPPOCAMPUS_l.nii.gz 3 1'%first_4d)
            os.system('fslroi %s PALLIDUM_l.nii.gz 4 1'%first_4d)
            os.system('fslroi %s PUTAMEN_l.nii.gz 5 1'%first_4d)
            os.system('fslroi %s THALAMUS_l.nii.gz 6 1'%first_4d)

            os.system('fslroi %s NACC_r.nii.gz 7 1'%first_4d)
            os.system('fslroi %s AMYGDALA_r.nii.gz 8 1'%first_4d)
            os.system('fslroi %s CAUDATE_r.nii.gz 9 1'%first_4d)
            os.system('fslroi %s HIPPOCAMPUS_r.nii.gz 10 1'%first_4d)
            os.system('fslroi %s PALLIDUM_r.nii.gz 11 1'%first_4d)
            os.system('fslroi %s PUTAMEN_r.nii.gz 12 1'%first_4d)
            os.system('fslroi %s THALAMUS_r.nii.gz 13 1'%first_4d)
            os.system('fslroi %s MIDBRAIN.nii.gz 14 1'%first_4d)

            os.system('fslmaths NACC_l -add CAUDATE_l -add PALLIDUM_l -add PUTAMEN_l -add THALAMUS_l -add THALAMUS_l -bin left_subcortical')
            os.system('fslmaths NACC_r -add CAUDATE_r -add PALLIDUM_r -add PUTAMEN_r -add THALAMUS_r -add THALAMUS_r -bin right_subcortical')
            os.system('fslmaths left_subcortical -add right_subcortical -bin subcortical')
            os.system('flirt -in subcortical.nii.gz -ref %s -out SUBCORTICAL_FUNC.nii.gz -applyxfm -init ANAT2FUNC_XFM.mat'%func_pproc_mean)
            os.system('fslmaths SUBCORTICAL_FUNC.nii.gz -thr 0.3 -bin ../NATIVE_FUNC_FIRST.nii.gz')

            # os.system('flirt -in PUTAMEN_l.nii.gz -ref %s -out func_PUTAMEN_l.nii.gz -applyxfm -init ANAT2FUNC_XFM.mat'%func_pproc_mean,)
            # os.system('flirt -in CAUDATE_l.nii.gz -ref %s -out func_CAUDATE_l.nii.gz -applyxfm -init ANAT2FUNC_XFM.mat'%func_pproc_mean,)
            # os.system('flirt -in PALLIDUM_l.nii.gz -ref %s -out func_PALLIDUM_l.nii.gz -applyxfm -init ANAT2FUNC_XFM.mat'%func_pproc_mean,)
            # os.system('flirt -in THALAMUS_l.nii.gz -ref %s -out func_THALAMUS_l.nii.gz -applyxfm -init ANAT2FUNC_XFM.mat'%func_pproc_mean,)
            #
            # os.system('fslmaths func_PUTAMEN_l.nii.gz -thr 0.3 -fillh -bin ../NATIVE_FUNC_PUTAMEN_L.nii.gz')
            # os.system('fslmaths func_CAUDATE_l.nii.gz -thr 0.3 -fillh -bin ../NATIVE_FUNC_CAUDATE_L.nii.gz')
            # os.system('fslmaths func_PALLIDUM_l.nii.gz -thr 0.3 -fillh -bin ../NATIVE_FUNC_PALLIDUM_L.nii.gz')

        ############################################################################################################
        # ANAT MASKS to MNI
        os.chdir(func2mni_dir)
        print '.....Transforming anatomical masks to MNI_2mm space'
        if not os.path.isfile('../MNI2mm_FUNC_GM.nii.gz'):
            os.system('WarpImageMultiTransform 3 %s MNI2mm_GM.nii.gz -R %s %s  %s'%(gm, mni_brain_2mm, warp, affine))
            os.system('WarpImageMultiTransform 3 %s MNI2mm_WM.nii.gz -R %s %s %s'%(wm, mni_brain_2mm, warp, affine))
            os.system('WarpImageMultiTransform 3 %s MNI2mm_CSF.nii.gz -R %s %s %s'%(csf, mni_brain_2mm, warp, affine))

            os.system('fslmaths MNI2mm_GM.nii.gz -mul ../REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz  -bin ../MNI2mm_FUNC_GM.nii.gz')
            os.system('fslmaths MNI2mm_WM.nii.gz -thr 0.99 -mul ../REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz  -bin ../MNI2mm_FUNC_WM.nii.gz')
            os.system('fslmaths MNI2mm_CSF.nii.gz -thr 0.99 -mul ../REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz  -bin ../MNI2mm_FUNC_CSF.nii.gz')

        if not os.path.isfile('../MNI2mm_FUNC_FIRST.nii.gz'):
            os.system('WarpImageMultiTransform 3 ../FUNC2ANAT/subcortical.nii.gz MNI2mm_SUBCORTICAL.nii.gz -R %s %s %s'%(mni_brain_2mm, warp, affine))
            os.system('fslmaths MNI2mm_SUBCORTICAL.nii.gz -mul ../REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz  -bin ../MNI2mm_FUNC_FIRST.nii.gz')

        if not os.path.isfile('../MNI2mm_FUNC_THALAMUS_L.nii.gz'):
            os.system('WarpImageMultiTransform 3 ../FUNC2ANAT/PUTAMEN_l.nii.gz MNI2mm_PUTAMEN_L.nii.gz -R %s %s %s'%(mni_brain_2mm, warp, affine))
            os.system('fslmaths MNI2mm_PUTAMEN_L.nii.gz -mul ../REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz  -fillh -bin ../MNI2mm_FUNC_PUTAMEN_L.nii.gz')

            os.system('WarpImageMultiTransform 3 ../FUNC2ANAT/CAUDATE_l.nii.gz MNI2mm_CAUDATE_L.nii.gz -R %s %s %s'%(mni_brain_2mm, warp, affine))
            os.system('fslmaths MNI2mm_CAUDATE_L.nii.gz -mul ../REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz  -fillh -bin ../MNI2mm_FUNC_CAUDATE_L.nii.gz')

            os.system('WarpImageMultiTransform 3 ../FUNC2ANAT/PALLIDUM_l.nii.gz MNI2mm_PALLIDUM_L.nii.gz -R %s %s %s'%(mni_brain_2mm, warp, affine))
            os.system('fslmaths MNI2mm_PALLIDUM_L.nii.gz -mul ../REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz  -fillh -bin ../MNI2mm_FUNC_PALLIDUM_L.nii.gz')

            os.system('WarpImageMultiTransform 3 ../FUNC2ANAT/THALAMUS_l.nii.gz MNI2mm_THALAMUS_L.nii.gz -R %s %s %s'%(mni_brain_2mm, warp, affine))
            os.system('fslmaths MNI2mm_THALAMUS_L.nii.gz -mul ../REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz  -fillh -bin ../MNI2mm_FUNC_THALAMUS_L.nii.gz')

        ############################################################################################################
        # FREESURFER MASKS to MNI
        freesurfer_mask_dir =  os.path.join(subject_dir, 'FUNC_TRANSFORM/FREESURFER_MASKS')
        mkdir_path(freesurfer_mask_dir)
        os.chdir(freesurfer_mask_dir)

        os.system('export SUBJECTS_DIR=%s'%(freesurfer_dir))

        mgz_t1  = os.path.join(freesurfer_dir, subject, 'mri', 'T1.mgz')
        mgz_seg = os.path.join(freesurfer_dir, subject, 'mri', 'aparc.a2009s+aseg.mgz')

        if not os.path.isfile('freesufer_2_native_space.mat'):
            #Convert freesurfer images to nifti rpi
            os.system('mri_convert %s freesurfer_aparc_aseg.nii.gz'%mgz_seg)
            os.system('mri_convert %s freesurfer_T1.nii.gz'%mgz_t1)
            os.system('fslswapdim freesurfer_aparc_aseg.nii.gz RL PA IS freesurfer_aparc_aseg_rpi.nii.gz')
            os.system('fslswapdim freesurfer_T1.nii.gz RL PA IS freesurfer_T1_rpi.nii.gz')

            # Register freesurfer T1/seg to MP2RAGE native space
            os.system('flirt -in freesurfer_T1_rpi.nii.gz -ref %s -omat freesufer_2_native_space.mat -dof 6 -out freesurfer_t1_native_space -cost mutualinfo' %anat)

        if not os.path.isfile('../MNI2mm_FUNC_INSULA_LA.nii.gz'):
            # Extract seg labels and transform to MNI space
            os.system('fslmaths freesurfer_aparc_aseg_rpi.nii.gz -thr 11148 -uthr 11148 -dilM -bin freesurfer_insula_anterior_left_rpi.nii.gz')
            os.system('flirt -in freesurfer_insula_anterior_left_rpi.nii.gz -ref %s -init freesufer_2_native_space.mat -applyxfm -out freesurfer_insula_anterior_left_rpi_native_space.nii.gz' %anat)

            #transform to MNI space
            os.system('WarpImageMultiTransform 3 freesurfer_insula_anterior_left_rpi_native_space.nii.gz MNI2mm_INSULA_anterior_left.nii.gz -R %s %s %s'%(mni_brain_2mm, warp, affine))
            os.system('fslmaths MNI2mm_INSULA_anterior_left.nii.gz -thr 0.5 -bin ../MNI2mm_FUNC_INSULA_LA.nii.gz')



run_functional_transform(['HR8T'], controls_datadir_a, workspace_a, freesurfer_dir_a)
#run_functional_transform(controls_a, controls_datadir_a, workspace_a, freesurfer_dir_a)
#run_functional_transform(patients_a, patients_datadir_a, workspace_a, freesurfer_dir_a)
#run_functional_transform(patients_b, patients_datadir_b, workspace_b, freesurfer_dir_b)
