__author__ = 'kanaan 17.12.2015'

import os
from variables.subject_list import *
from utilities.utils import mkdir_path
from utilities.nuisance import *
from utilities.bandpass import *
import shutil
import subprocess
import commands



def run_randomise_baseline(population_1, population_2, workspace_dir):

    print '#######################################'

    for string in ecm_strings:
        print 'Running FSL-randomise on  ECM maps = %s'%string[16:]
        print ''
        ecm_list_1 = []
        ecm_list_2 = []
        subs_1     = []
        subs_2     = []

        for subject in population_1:
            ecm = os.path.join(workspace_dir, subject, 'FAST_ECM/%s/zscore_fastECM.nii.gz'%string)
            if os.path.isfile(ecm):
                ecm_list_1.append(ecm)
                subs_1.append(subject)

        for subject in population_2:
            ecm = os.path.join(workspace_dir, subject, 'FAST_ECM/%s/zscore_fastECM.nii.gz'%string)
            if os.path.isfile(ecm):
                ecm_list_2.append(ecm)
                subs_2.append(subject)

        print '...Concatenating ECM maps into 4D'
        # calculate group mean
        out_dir = os.path.join(workspace_dir, 'STATISTICS/ECM/%s'%string)
        mkdir_path(out_dir)
        os.chdir(out_dir)
        if not os.path.isfile('ECM_all_concat.nii.gz'):
            os.system('fslmerge -t ECM_controls_concat.nii.gz %s'%' '.join(ecm_list_1))
            print 'Controls n=%s'%len(subs_1), subs_1
            os.system('fslmerge -t ECM_patients_concat.nii.gz %s'%' '.join(ecm_list_2))
            print 'Patients n=%s'%len(subs_2), subs_2
            print ''
            os.system('fslmaths ECM_controls_concat.nii.gz -Tmean ECM_controls_mean.nii.gz')
            os.system('fslmaths ECM_patients_concat.nii.gz -Tmean ECM_patients_mean.nii.gz')
            os.system('fslmerge -t ECM_all_concat.nii.gz ECM_controls_concat.nii.gz ECM_patients_concat.nii.gz')

        print '...Running Non-paramteric permuation tests'
        # Runs FSL Randomise - nonparametric permutation inference
        # Two-Sample Unpaired T-test

        group_gm_mask     = '/SCR4/workspace/project_GluRest/OUT_DIR_A/GluConnectivity/GM2mm_bin.nii'

        glm_mat = os.path.join(workspace_dir, 'STATISTICS/ECM/GLM_baseline/randomise_baseline.mat')
        glm_con = os.path.join(workspace_dir, 'STATISTICS/ECM/GLM_baseline/randomise_baseline.con')
        if not os.path.isfile('randomise_baseline_tfce_corrp_tstat2.nii.gz'):
            os.system('randomise -i ECM_all_concat.nii.gz -o randomise_baseline -d %s -t %s -T -R -P -N --uncorrp -n 5000 -m %s '
                      %(glm_mat, glm_con, group_gm_mask))

run_randomise_baseline(controls_a_qc, patients_a_qc, workspace_a)

