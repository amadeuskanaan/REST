__author__ = 'kanaan 03.12.2015'

import os
from variables.subject_list import *
from utilities.utils import *
import shutil
import nipype.interfaces.ants as ants


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
                                       Functional Pre-processing
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def run_functional_preprocessing(population, datadir, workspace_dir):


     for subject in population:
        print '###############################################################################'
        print 'Processing Functional Data for subject %s' %subject
        print ''

        # input
        rest        =  os.path.join(datadir, subject, 'NIFTI', 'REST.nii')
        rest_rl     =  os.path.join(datadir, subject, 'NIFTI', 'REST_SE.nii')
        rest_lr     =  os.path.join(datadir, subject, 'NIFTI', 'REST_SE_INVPOL.nii')

        #output
        func_dir =  os.path.join(workspace_dir, subject, 'FUNC_PPROC')
        mkdir_path(func_dir)
        os.chdir(func_dir)


        ############################################################################################################
        # Reorient and Drop TRs

        print '.....Reorienting Data'
        if not os.path.isfile(os.path.join(func_dir, 'REST_DROP_RPI.nii.gz' )):
            # Deoblique data
            shutil.copy(rest, os.path.join(func_dir, 'REST.nii'))
            shutil.copy(rest_rl, os.path.join(func_dir, 'REST_SE.nii'))
            shutil.copy(rest_lr, os.path.join(func_dir, 'REST_SE_INVPOL.nii'))
            os.system('3drefit -deoblique REST.nii')
            os.system('3drefit -deoblique REST_SE.nii')
            os.system('3drefit -deoblique REST_SE_INVPOL.nii')

            # Drop 5TRs
	    if subject is 'RMNT':
	      f = '[5..420]'  
	    else:
	      f = '[5..421]'
            os.system('3dcalc -a REST.nii%s -expr "a" -prefix REST_DROP.nii.gz'%f)
		   
            # Reorient to RPI
            os.system('fslswapdim REST_DROP.nii.gz RL PA IS REST_DROP_RPI.nii.gz')
            os.system('bet REST_DROP_RPI.nii.gz REST_DROP_RPI_BRAIN.nii.gz -f 0.50 -F -m -t -g 0.00')
            os.system('fslswapdim REST_SE.nii RL PA IS REST_SE_RPI.nii')
            os.system('fslswapdim REST_SE_INVPOL.nii RL PA IS REST_SE_INVPOL_RPI.nii')

        ############################################################################################################
        # CALC DISTORTION CORRECTION
        print '.....Running Distortion Correction'

        disco_dir =  os.path.join(func_dir, 'DISTORTION_CORRECTION')
        mkdir_path(disco_dir)
        os.chdir(disco_dir)

        if not os.path.isfile(os.path.join(func_dir, 'REST_DISCO.nii.gz')):

            #1. merge blips
            os.system('fslmerge -t BLIPS_MERGED.nii.gz ../REST_SE_RPI.nii ../REST_SE_INVPOL_RPI.nii')

            #2. Calc topup field
            os.system(' '.join([ 'topup',
                                 '--config=b02b0.cnf',
                                 '--datain=/scr/sambesi1/workspace/Projects/GluREST/functional/datain.txt',
                                 '--imain=BLIPS_MERGED.nii.gz',
                                 '--out=TOPUP',
                                 '--iout=TOPUP_CORRECTED.nii.gz',
                                 '--fout=TOPUP_FIELD.nii.gz',
                                 '--logout=TOPUP.log']))

            #3. Apply topup
            os.system(' '.join([ 'applytopup',
                                 '--datain=/scr/sambesi1/workspace/Projects/GluREST/functional/datain.txt',
                                 '--imain=../REST_DROP_RPI.nii.gz',
                                 '--inindex=1',
                                 '--topup=TOPUP',
                                 '--out=REST_DISCO.nii.gz',
                                 '--method=jac']))

            os.system('cp REST_DISCO.nii.gz ../REST_DISCO.nii.gz')
        ############################################################################################################
        # CALC MOTION CORRECTION
        print '.....Running Motion Correction'

        moco_dir =  os.path.join(func_dir, 'MOTION_CORRECTION')
        mkdir_path(moco_dir)
        os.chdir(moco_dir)

        # 1. MOCO Step 1
        if not os.path.isfile(os.path.join(moco_dir, 'REST_DISCO_mean.nii.gz')):
            os.system('fslmaths ../REST_DISCO.nii.gz -Tmean REST_DISCO_mean.nii.gz')

        # AFNI
        if not os.path.isfile('./REST_DISCO_MOCO_1.nii.gz'):
            os.system('3dvolreg -Fourier -twopass -1Dfile REST_DISCO_MOCO_1.1D -1Dmatrix_save REST_DISCO_MOCO_1.aff12.1D '
                      '-prefix REST_DISCO_MOCO_1.nii.gz -base REST_DISCO_mean.nii.gz -zpad 4 -maxdisp1D REST_DISCO_MOCO_1_MX.1D '
                      '../REST_DISCO.nii.gz')
            os.system('fslmaths REST_DISCO_MOCO_1.nii.gz -Tmean REST_DISCO_MOCO_1_mean.nii.gz')


        if not os.path.isfile('./REST_DISCO_MOCO_2.nii.gz'):
            os.system('3dvolreg -Fourier -twopass -1Dfile REST_DISCO_MOCO_2.1D -1Dmatrix_save REST_DISCO_MOCO_2.aff12.1D '
                      '-prefix REST_DISCO_MOCO_2.nii.gz -base  REST_DISCO_MOCO_1_mean.nii.gz -zpad 4 -maxdisp1D REST_DICSO_MOCO_2_MX.1D '
                      '../REST_DISCO.nii.gz')
            os.system('fslmaths REST_DISCO_MOCO_2.nii.gz -Tmean REST_DISCO_MOCO_2_mean.nii.gz')

            mats = create_fsl_mats('REST_DISCO_MOCO_2.aff12.1D')

        # ############################################################################################################
        # ########deprecated for now
        # # APPLY UNIFORM WARP
        # print '.....Combining Warps and Affines in one Resampling Step'
        # #
        # unify_dir =  os.path.join(func_dir, 'CONVERT_XFM')
        # mkdir_path(unify_dir)
        # os.chdir(unify_dir)
        #
        # if not os.path.isfile('REST_UNIFED_WARP.nii.gz'):
        #     if not os.path.isfile('vol0416.nii.gz'):
        #         # split functional frames
        #         os.system('fslsplit ../REST_DROP_RPI.nii.gz -t moco_split')
        #
        #     if not os.path.isfile('WARP_FIELD.nii.gz'):
        #         # Multiply voxel shift may by epi readout time
        #         os.system('fslmaths  %s/TOPUP_FIELD.nii.gz -mul 0.0589587878 VOXEL_SHIFT_MAP'%disco_dir)
        #         os.system('convertwarp -r  %s/REST_DISCO_MOCO_2_mean.nii.gz -s VOXEL_SHIFT_MAP.nii.gz --relout --rel -d y -o WARP_FIELD.nii.gz'%moco_dir)
        #         # os.system('fugue --loadfmap=%s/TOPUP_FIELD.nii.gz --dwell=0.0589587878 --saveshift=SHIFT_MAP '%disco_dir)
        #         # os.system('convertwarp -r %s/REST_DISCO_MOCO_2_mean.nii.gz -s VOXEL_SHIFT_MAP.nii.gz --relout --rel -d y -o WARP_FIELD.nii.gz'%moco_dir)
        #
        #         #fugue --loadfmap=${WD}/FieldMap2${TXwImageBasename} --dwell=${TXwSampleSpacing} --saveshift=${WD}/FieldMap2${TXwImageBasename}_ShiftMap.nii.gz ^C
        #
        #     # APPLY WARP with Motion correction as premet
        #     for i in xrange(0,416):
        #         frame = '{0:0>4}'.format(i)
        #         os.system('applywarp -i vol%s -o applywarp%s.nii.gz -r ../MOTION_CORRECTION/REST_DISCO_MOCO_2_mean.nii.gz '
        #                   '--premat=../MOTION_CORRECTION/MATS/MAT_%s --rel -w WARP_FIELD.nii.gz --interp=spline -d float' %(frame,frame, frame))
        #         os.system('fslmerge -t REST_UNIFED_WARP.nii.gz applywarp*')

        ###########################################################################################################
        #Intensity Normalization
        print '.....Normalizing intensity to Mode 1000 and deskulling'
        os.chdir(func_dir)#
        if not os.path.isfile('REST_PPROC_NATIVE_BRAIN.nii.gz'):
            os.system('fslmaths MOTION_CORRECTION/REST_DISCO_MOCO_2.nii.gz -ing 1000 REST_PPROC_NATIVE.nii.gz')
            os.system('bet REST_PPROC_NATIVE.nii.gz REST_PPROC_NATIVE_BRAIN.nii.gz -f 0.50 -F -m -t -g 0.00')
            os.system('3dTstat -mean -prefix REST_PPROC_NATIVE_BRAIN_mean.nii.gz  REST_PPROC_NATIVE_BRAIN.nii.gz ')
            os.system('fslmaths REST_PPROC_NATIVE_BRAIN_mask.nii.gz -ero -ero REST_PPROC_NATIVE_BRAIN_mask_ero.nii.gz')

run_functional_preprocessing(['HCTT'], controls_datadir_a, workspace_a)
#run_functional_preprocessing(controls_a, controls_datadir_a, workspace_a)
#run_functional_preprocessing(patients_a, patients_datadir_a, workspace_a)
#run_functional_preprocessing(patients_b, patients_datadir_b, workspace_b)
