__author__ = 'kanaan 05.01.2016'

import os
from utilities.utils import *
from variables.subject_list import *


def run_randomise_sca_baseline(population_1, population_2, workspace_dir, mask_name):

    for string in resid_strings:
        sca_list_1 = []
        sca_list_2 = []
        subs_1     = []
        subs_2     = []

        for subject in population_1:
            sca = os.path.join(workspace_dir, subject, 'SCA/%s/SCA_Z_%s.nii.gz'%(mask_name, string))
            if os.path.isfile(sca):
                sca_list_1.append(sca)
                subs_1.append(subject)

        for subject in population_2:
            sca = os.path.join(workspace_dir, subject, 'SCA/%s/SCA_Z_%s.nii.gz'%(mask_name, string))
            if os.path.isfile(sca):
                sca_list_2.append(sca)
                subs_2.append(subject)

        print '...Concatenating SCA maps into 4D'
        out_dir = os.path.join(workspace_dir, 'STATISTICS/SCA/%s/%s'%(mask_name, string))
        mkdir_path(out_dir)
        os.chdir(out_dir)

        if not os.path.isfile('SCA_all_concat.nii.gz'):
            os.system('fslmerge -t SCA_controls_concat.nii.gz %s'%' '.join(sca_list_1))
            print 'Controls n=%s'%len(subs_1), subs_1

            os.system('fslmerge -t SCA_patients_concat.nii.gz %s'%' '.join(sca_list_2))
            print 'Patients n=%s'%len(subs_2), subs_2
            print ''

            os.system('fslmaths SCA_controls_concat.nii.gz -Tmean SCA_controls_mean.nii.gz')
            os.system('fslmaths SCA_patients_concat.nii.gz -Tmean SCA_patients_mean.nii.gz')
            os.system('fslmerge -t SCA_all_concat.nii.gz SCA_controls_concat.nii.gz SCA_patients_concat.nii.gz')

        print '...Running Non-paramteric permuation tests'
        # Runs FSL Randomise - nonparametric permutation inference
        # Two-Sample Unpaired T-test
        glm_mat = os.path.join(workspace_dir, 'STATISTICS/ECM/GLM_baseline/randomise_baseline.mat')
        glm_con = os.path.join(workspace_dir, 'STATISTICS/ECM/GLM_baseline/randomise_baseline.con')

        if not os.path.isfile('randomise_baseline_tfce_corrp_tstat2.nii.gz'):
            os.system('randomise -i SCA_all_concat.nii.gz -o randomise_baseline -d %s -t %s -T -R -P -N --uncorrp -n 500'
                      %(glm_mat, glm_con))

# run_randomise_sca_baseline(controls_a_qc, patients_a_qc, workspace_a, 'PUTAMEN')
run_randomise_sca_baseline(controls_a_qc, patients_a_qc, workspace_a, 'CAUDATE')
# run_randomise_sca_baseline(controls_a_qc, patients_a_qc, workspace_a, 'PALLIDUM')
# run_randomise_sca_baseline(controls_a_qc, patients_a_qc, workspace_a, 'THALAMUS')
# run_randomise_sca_baseline(controls_a_qc, patients_a_qc, workspace_a, 'FIRST')
