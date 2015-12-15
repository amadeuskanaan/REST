__author__ = 'kanaan'

import os
import shutil
from variables.subject_list import *
from utilities.utils import mkdir_path

def grab_subcortical_optimized_tissue_masks(population, data_dir, population_name, data_dumpdir):

    print '#############################################################################'
    print ''
    print '                 RUNNNING PROJECT NMR-093%s %s' %(data_dumpdir[-10], population_name)
    print ''
    print '#############################################################################'

    count = 0
    for subject in population:
        count +=1
        print '%s-Copying preprocessed anatomical data for %s %s-%s'%(count,population_name[:-1], subject, data_dumpdir[-10])

        #inputs
        data_in_dir = os.path.join(data_dir, population_name, subject, )
        gm = os.path.join(data_in_dir,'segmentation_spm', 'TISSUE_CLASS_1_GM_OPTIMIZED.nii.gz')
        wm = os.path.join(data_in_dir,'segmentation_spm', 'TISSUE_CLASS_2_WM_OPTIMIZED.nii.gz')
        cm = os.path.join(data_in_dir,'segmentation_spm', 'TISSUE_CLASS_3_CSF_OPTIMIZED.nii.gz')
        anat = os.path.join(data_in_dir,'anatomical_original', 'ANATOMICAL_DESKULL_RPI.nii.gz')
        first = os.path.join(data_in_dir,'segmentation_spm', 'FIRST_subcortical', 'FIRST_all_fast_firstseg.nii.gz')
        firsto = os.path.join(data_in_dir,'segmentation_spm', 'FIRST_subcortical', 'FIRST_all_fast_origsegs.nii.gz')

        #outputs
        #mkdir_path(os.path.join(data_dumpdir, subject, 'NIFTI'))
        data_out_dir = os.path.join(data_dumpdir, subject, 'NIFTI')

        # copy segments to out dir
        shutil.copy(gm, data_out_dir)
        shutil.copy(wm, data_out_dir)
        shutil.copy(cm, data_out_dir)
        shutil.copy(anat, os.path.join(data_out_dir, 'MP2RAGE_DESKULL_RPI.nii.gz'))
        shutil.copy(first, os.path.join(data_out_dir, 'FIRST_3d.nii.gz'))
        shutil.copy(firsto, os.path.join(data_out_dir, 'FIRST_4d.nii.gz'))

if __name__ == "__main__":
    # grab_subcortical_optimized_tissue_masks(population= ['GHAT'], data_dir=mrs_datadir_a, population_name= 'controls', data_dumpdir= controls_datadir_a)
    grab_subcortical_optimized_tissue_masks(population= controls_a, data_dir=mrs_datadir_a, population_name= 'controls', data_dumpdir= controls_datadir_a)
    grab_subcortical_optimized_tissue_masks(population= patients_a, data_dir=mrs_datadir_a, population_name= 'patients', data_dumpdir= patients_datadir_a)
    grab_subcortical_optimized_tissue_masks(population= controls_b, data_dir=mrs_datadir_b, population_name= 'controls', data_dumpdir= controls_datadir_b)
    grab_subcortical_optimized_tissue_masks(population= patients_b, data_dir=mrs_datadir_b, population_name= 'patients', data_dumpdir= patients_datadir_b)
