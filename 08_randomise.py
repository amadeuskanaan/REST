__author__ = 'kanaan 17.12.2015'

import os
from variables.subject_list import *
from utilities.utils import mkdir_path
from utilities.nuisance import *
from utilities.bandpass import *
import shutil
import subprocess
import commands

def create_ecm_group_mean(population, population_name, workspace_dir):

    ecm_list = []

    ecm_strings = ["RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm",
                  "RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm",
                  "RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm",
                  "RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp",
                  "RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp",
                  "RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp"]
    print '#######################################'
    print 'Creating output files for FSL-randomise'
    print ''
    print 'Grabbing Files '
    for string in ecm_strings:
        for subject in population:
            subject_dir    = os.path.join(workspace_dir, subject)
            ecm            = os.path.join(subject_dir, 'FAST_ECM/%s/zscore_fastECM.nii.gz'%string)
            if os.path.isfile(ecm):
                ecm_list.append(ecm)

        print '-Concatenating maps and calculating mean for %s --->%s'%(population_name, string[16:])
        # calculate group mean
        out_dir = os.path.join(workspace_dir, 'STATISTICS/ECM/%s'%string)
        mkdir_path(out_dir)
        os.chdir(out_dir)
        os.system('fslmerge -t ECM_%s_concat.nii.gz %s'%(population_name, ' '.join(ecm_list)))
        os.system('fslmaths ECM_%s_concat.nii.gz -Tmean ECM_%s_mean.nii.gz'%(population_name, population_name))


create_ecm_group_mean(controls_a, 'controls', workspace_a)
create_ecm_group_mean(patients_a, 'patients', workspace_a)



# viz
from surfer import Brain, io

# define stat file and freesurfer reg file
ecm_file = '/scr/sambesi4/workspace/project_REST/study_a/STATISTICS/ECM/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp/ECM_controls_mean.nii.gz'
reg_file = '/afs/cbs.mpg.de/software/freesurfer/5.3.0/ubuntu-precise-amd64/average/mni152.register.dat'

#project data to fsaverage
surf_data = io.project_volume_data(ecm_file, "lh", reg_file)

brain = Brain("fsaverage", "lh", "pial")
#brain.add_data(surf_data* 10, min = -24, max = 24 , colormap="jet")
brain.add_data(surf_data *10 ,min = -5, max = 5,   colormap="jet")



#
# # viz
# from surfer import Brain, io
#
# # define stat file and freesurfer reg file
# ecm_file_c = '/scr/sambesi4/workspace/project_REST/study_a/STATISTICS/ECM/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp/ECM_controls_mean.nii.gz'
# ecm_file_p = '/scr/sambesi4/workspace/project_REST/study_a/STATISTICS/ECM/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp/ECM_patients_mean.nii.gz'
# ecm_file_diff = '/scr/sambesi4/workspace/project_REST/study_a/STATISTICS/ECM/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp/pminusc.nii.gz'
# reg_file = '/afs/cbs.mpg.de/software/freesurfer/5.3.0/ubuntu-precise-amd64/average/mni152.register.dat'
#
# #project data to fsaverage
# surf_data = io.project_volume_data(ecm_file_diff, "lh", reg_file)
#
# brain = Brain("fsaverage", "lh", "pial")
# #brain.add_data(surf_data* 10, min = -24, max = 24 , colormap="jet")
# brain.add_data(surf_data *10 ,min = -5, max = 5,   colormap="jet")

# #####
# # define stat file and freesurfer reg file
# from surfer import Brain, io
# ecm_file = '/scr/sambesi4/workspace/project_REST/study_a/STATISTICS/ECM/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp/ECM_controls_mean.nii.gz'
# reg_file = '/afs/cbs.mpg.de/software/freesurfer/5.3.0/ubuntu-precise-amd64/average/mni152.register.dat'
# brain = Brain("fsaverage", "split", "pial",  size=(800, 400))
# surf_l = io.project_volume_data(ecm_file, 'lh', reg_file)
# surf_r = io.project_volume_data(ecm_file, 'rh', reg_file)
#
# brain.add_data(surf_l*10,min=-5, max = 5, colormap='jet', hemi = 'lh')
# brain.add_data(surf_r*10,min=-5, max = 5, colormap='jet', hemi = 'rh')
