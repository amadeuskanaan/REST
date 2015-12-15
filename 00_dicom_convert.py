__author__ = 'kanaan 08.03.2015'

import os
import sys
import shutil
import subprocess
from nipype.interfaces.mipav.developer import JistIntensityMp2rageMasking
from nipype.interfaces.mipav.developer import MedicAlgorithmSPECTRE2010
import nipype.interfaces.fsl as fsl
from variables.subject_list import *
from utilities.utils import locate

def convert_scanner_data(population, afs_dir, data_dumpdir):

    count=0
    for subject in population:
        count +=1
        print '========================================================================================'
        print      '%s- Dicom Conversion and Anatomical Preprocessing for subject %s_%s' %(count,subject, afs_dir[-2])
        print '========================================================================================'

        for folder in os.listdir(afs_dir):
            if folder.startswith('p'):

                '===================================================================================================='
                '                                        DICOM to NIFTI                                              '
                '===================================================================================================='

                #set dicom dir
                if os.path.isdir(os.path.join(afs_dir, folder, subject, 'DICOM')):
                    dicom_dir =  os.path.join(afs_dir, folder, subject, 'DICOM')

                    if os.path.isdir(os.path.join(data_dumpdir, subject)):
                         pass
                    else:
                         os.makedirs(os.path.join(data_dumpdir, subject))
                    subject_dir = os.path.join(data_dumpdir, subject)


                    # create output out dir
                    try:
                        os.makedirs(os.path.join(subject_dir, 'NIFTI'))
                    except OSError:
                        nifti_dir  = str(os.path.join(subject_dir, 'NIFTI'))
                    nifti_dir      = str(os.path.join(subject_dir, 'NIFTI'))

                    # convert dicoms to niftis
                    # ensure conversion hasnt been completed before
                    if os.path.isfile(os.path.join(nifti_dir, 'REST.nii')):
                        print 'Dicom Conversion already completed...... moving on'
                    else:
                        print 'Converting DICOM to NIFTI'
                        convert_cmd = ['isisconv',
                                       '-in',
                                       '%s' %dicom_dir,
                                       '-out',
                                       '%s/%s_S{sequenceNumber}_{sequenceDescription}_{echoTime}.nii' %(nifti_dir, subject),
                                       '-rf',
                                       'dcm',
                                       '-wdialect',
                                       'fsl']
                        print subprocess.list2cmdline(convert_cmd)
                        subprocess.call(convert_cmd)


                        #rename outputs

                        for file in os.listdir(nifti_dir):
                            if 'mp2rage_p3_602B_INV1_2.98' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'MP2RAGE_INV1.nii')))
                            elif 'mp2rage_p3_602B_INV2_2.98' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'MP2RAGE_INV2.nii')))
                            elif 'mp2rage_p3_602B_DIV_Images_2.98' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'MP2RAGE_DIV.nii')))
                            elif 'mp2rage_p3_602B_T1_Images_2.98' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'MP2RAGE_T1MAPS.nii')))
                            elif 'mp2rage_p3_602B_UNI_Images_2.98' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'MP2RAGE_UNI.nii')))
                            elif 'resting' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'REST.nii')))
                            elif 'mbep2d_se_52' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'REST_SE.nii')))
                            elif 'se_invpol_52' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'REST_SE_INVPOL.nii')))
                            elif 'bvec' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'DWI_BVEC.bvec')))
                            elif 'AP_unwarp_diff' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'DWI_AP.nii')))
                            elif 'PA_unwarp_diff' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'DWI_PA.nii')))
                            elif 'mbep2d_diff_80' in file:
                                os.rename(str(os.path.join(nifti_dir, file)),
                                          str(os.path.join(nifti_dir, 'DWI.nii')))

                            # remove irrelevent files to conserve space
                            irrelvent_lst = ['AAH', 'AX', 'ax', 'Cor', 'cor', 'COR', 'hip',
                                             'Hip', 'slab', 'Modus', 'acpc', 'DUMMY', 'dummy',
                                             'short', 'SLAB' ]
                            try:
                                for string in irrelvent_lst:
                                    if string in file:
                                        os.remove(str(os.path.join(nifti_dir, file)))
                            except OSError:
                                print 'cant delete file %s' %str(os.path.join(nifti_dir, file))

                    '===================================================================================================='
                    '                                  Denoising MPRAGE Anatomical                                       '
                    '===================================================================================================='

                    if os.path.isfile(os.path.join(nifti_dir, 'MP2RAGE_BRAIN.nii')):
                        print 'MP2RAGE already deskulled............... moving on'
                    else:
                        print 'Deskulling  MP2RAGE'

                        try:
                            mp2rage_uni        = locate('MP2RAGE_UNI.nii', nifti_dir)
                            mp2rage_inv2       = locate('MP2RAGE_INV2.nii', nifti_dir)
                            mp2rage_t1maps     = locate('MP2RAGE_T1MAPS.nii', nifti_dir)
                        except OSError:
                            continue

                        try:
                            os.makedirs(os.path.join(nifti_dir, 'MIPAV_OUTPUTS'))
                        except OSError:
                            mipav_dir  = str(os.path.join(nifti_dir, 'MIPAV_OUTPUTS'))
                        mipav_dir  = str(os.path.join(nifti_dir, 'MIPAV_OUTPUTS'))

                        os.chdir(mipav_dir)

                        t1_threshold                 = fsl.Threshold()
                        t1_threshold.inputs.in_file  = mp2rage_t1maps
                        t1_threshold.inputs.thresh   = 1
                        t1_threshold.inputs.args     = '-uthr 4000'
                        t1_threshold.run()

                        mp2rage_t1maps_thr = locate('MP2RAGE_T1MAPS_thresh.nii.gz', mipav_dir)

                        anat_denoise  = JistIntensityMp2rageMasking(outMasked=True,outMasked2=True,outSignal2=True,inSkip = 'true')
                        anat_denoise.inputs.inQuantitative = mp2rage_t1maps_thr
                        anat_denoise.inputs.inT1weighted   = mp2rage_uni
                        anat_denoise.inputs.inSecond       = mp2rage_inv2
                        anat_denoise.run()

                        uni_denoised   = locate('outMasked2.nii', mipav_dir)

                        '===================================================================================================='
                        '                                          Deskulling                                                '
                        '===================================================================================================='

                        anat_deskull   = MedicAlgorithmSPECTRE2010(inAtlas       = str('/afs/cbs.mpg.de/software/cbstools/3.0/jist-cruise/Atlas/spectre/oasis-3-v2.txt'),
                                                                   #inInitial    = 5,
                                                                   #inInitial2   = 0.35 ,
                                                                   #inMinimum    = 0.1,
                                                                   #inSmoothing  = 0.02,
                                                                   #inBackground = 0.001,
                                                                   outOriginal   = True,
                                                                   #inOutput     = 'true',
                                                                   outStripped   = True,
                                                                   outMask       = True,
                                                                   inFind        = 'true',
                                                                   xMaxProcess   = 0,
                                                                   inMMC         = 2,
                                                                   inMMC2        = 2)
                        anat_deskull.inputs.inInput = uni_denoised
                        anat_deskull.run()

                        shutil.move(str(os.path.join(mipav_dir, 'outStripped.nii')),
                                    str(os.path.join(nifti_dir, 'MP2RAGE_BRAIN.nii')))

                    if os.path.isfile(os.path.join(nifti_dir, 'MP2RAGE_BRAIN.nii')):
                        brain = os.path.join(nifti_dir, 'MP2RAGE_BRAIN.nii')
                        print 'Path =  %s' %brain

    print '========================================================================================'

if __name__ == "__main__":
    convert_scanner_data(controls_a, afsdir_a, controls_datadir_a)
    convert_scanner_data(patients_a, afsdir_a, patients_datadir_a)
    convert_scanner_data(controls_b, afsdir_b, controls_datadir_b)
    convert_scanner_data(patients_b, afsdir_b, patients_datadir_b)
