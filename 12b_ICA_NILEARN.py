__author__ = 'kanaan 15.02.2016'

import os
from utilities.utils import mkdir_path
from variables.subject_list import *
from nipype.interfaces.fsl import Merge
from nipype.interfaces import afni as afni
import nibabel as nb
import commands
from nilearn.decomposition.canica import CanICA


def run_group_ica(population, workspace_dir):

    print '#############################################################################'
    print ''
    print '                          RUNNNING GROUP ICA'
    print ''
    print '#############################################################################'

    #define ica outputdir

    group_ica_dir = os.path.join(workspace_dir, 'NILEARN_GROUP_ICA')
    mkdir_path(group_ica_dir)
    os.chdir(group_ica_dir)

    #concatenate fmri preprocessed data.
    preprocssed_all =[]
    for subject in population:
        subject_dir = os.path.join(workspace_dir, subject)
        preprocssed_2mm = os.path.join(subject_dir, 'FUNC_TRANSFORM/ICAready_REST_PPROC_MNI2mm_fwhm_bp.nii.gz')
        preprocssed_4mm = os.path.join(subject_dir, 'FUNC_TRANSFORM/ICAready_REST_PPROC_MNI2mm_fwhm_bp_4mm.nii.gz')
        if not os.path.isfile(preprocssed_4mm):
            print 'resampleing subject %s to 4mm'%subject
            os.system('flirt -in %s -ref %s -out %s/FUNC_TRANSFORM/ICAready_REST_PPROC_MNI2mm_fwhm_bp_4mm.nii.gz -applyisoxfm 4'%(preprocssed_2mm, mni_brain_2mm, subject_dir))
            preprocssed_4mm = os.path.join(subject_dir, 'FUNC_TRANSFORM/ICAready_REST_PPROC_MNI2mm_fwhm_bp_4mm.nii.gz')
        print preprocssed_4mm
   
        if os.path.isfile(preprocssed_4mm):
            preprocssed_all.append(preprocssed_4mm)
        else:
            print 'subjects with missing data ',subject

    #preprocssed_all = ['/scr/sambesi4/workspace/project_REST/study_a/BM8X/FUNC_TRANSFORM/ICAready_REST_PPROC_MNI2mm_fwhm_bp_4mm.nii.gz']
    #print preprocssed_all


run_group_ica(population_qc, workspace_a)


