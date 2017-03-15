__author__ = 'kanaan 15.02.2016'

import os
from utilities.utils import mkdir_path
from variables.subject_list import *
from nipype.interfaces.fsl import Merge
from nipype.interfaces import afni as afni
import nibabel as nb
import commands



def run_group_ica(population, workspace_dir):

    print '#############################################################################'
    print ''
    print '                          RUNNNING GROUP ICA'
    print ''
    print '#############################################################################'

    #define ica outputdir

    group_ica_dir = os.path.join(workspace_dir, 'MELODIC_GROUP_ICA')
    mkdir_path(group_ica_dir)
    os.chdir(group_ica_dir)

    # concatenate fmri preprocessed data.
    preprocssed_all =[]
    for subject in population:
        subject_dir = os.path.join(workspace_dir, subject)
        preprocssed_2mm = os.path.join(subject_dir, 'FUNC_TRANSFORM/ICAready_REST_PPROC_MNI2mm_fwhm_bp.nii.gz')
        preprocssed_4mm = os.path.join(subject_dir, 'FUNC_TRANSFORM/ICAready_REST_PPROC_MNI2mm_fwhm_bp_4mm.nii.gz')
        if not os.path.isfile(preprocssed_4mm):
            print 'resampleing subject %s to 4mm'%subject
            os.system('flirt -in %s -ref %s -out %s/FUNC_TRANSFORM/ICAready_REST_PPROC_MNI2mm_fwhm_bp_4mm.nii.gz -applyisoxfm 4'%(preprocssed_2mm, mni_brain_2mm, subject_dir))
            preprocssed_4mm = os.path.join(subject_dir, 'FUNC_TRANSFORM/ICAready_REST_PPROC_MNI2mm_fwhm_bp_4mm.nii.gz')
        #print preprocssed_4mm

        if os.path.isfile(preprocssed_4mm):
            preprocssed_all.append(preprocssed_4mm)
        else:
            print 'subjects with missing data ',subject

    input_file = os.path.join(group_ica_dir, 'input_list.txt')
    with open(input_file, 'w') as file:
        for i in preprocssed_all:
            file.write(i+'\n')
    print input_file

    def run_melodic(func, brain_mask, TR, melodic_dir):

        # Run MELODIC
        #os.system('melodic --in=%s --outdir=%s --mask=%s -Ostats --nobet --mmthresh=0.5 --report --tr=%s' %(func, melodic_dir, brain_mask, str(TR)))
        os.system(' '.join([ 'melodic',
                             '--in=' + func,
                             '--mask=' + brain_mask,
                             '-v',
			                 '-d 25',
                             '--outdir='  + melodic_dir,
                             '--Ostats --nobet --mmthresh=0.5 --report',
                             '--tr=' + str(TR)]))
        if os.path.isfile(os.path.join(melodic_dir, 'melodic_IC.nii.gz')):
            # Get number of components
            melodic_4d = nb.load(os.path.join(melodic_dir, 'melodic_IC.nii.gz'))
            n_componenets = melodic_4d.shape[3]

            for n_componenet in range(1,n_componenets):
                z_thresh = os.path.join(melodic_dir, 'stats/thresh_zstat%s.nii.gz'%n_componenet)

                cmd = ' '.join(['fslinfo', z_thresh, '| grep dim4 | head -n1 | awk \'{print $2}\''])
                z_thresh_dim4 = int(float(commands.getoutput(cmd)))

                # Zero-pad the IC number and extract the 3D data........
                # For cases where the mixture modeling does not converge, 2nd img in the 4th dimension wil be the results of the null hypothesis test.
                cmd = ' '.join(['zeropad', str(n_componenet), '4'])
                z_thresh_zeropad = os.path.join(melodic_dir,'thr_zstat' + commands.getoutput(cmd))

                # Extract last spatial map within the thresh_zstat file
                os.system('fslroi %s %s %s 1' %(z_thresh, z_thresh_zeropad, str(z_thresh_dim4-1)))

            # Merge and subsequently remove all mixture modeled Z-maps within the output directory
            z_thresh_zeropad_all = os.path.join(melodic_dir, 'thr_zstat????.nii.gz')
            z_thresh_merged = os.path.join(melodic_dir, 'melodic_IC_thr.nii.gz')

            os.system('fslmerge -t %s %s ' %(z_thresh_merged , z_thresh_zeropad_all))
            os.system('rm ' + z_thresh_zeropad_all)

            # Apply the mask to the merged file (in case a melodic-directory was predefined and run with a different mask)
            os.system('fslmaths %s -mas %s %s' %(z_thresh_merged, brain_mask, z_thresh_merged))

    run_melodic(input_file, mni_brain_mask_4mm, TR = 1.4, melodic_dir = group_ica_dir)

run_group_ica(population_qc, workspace_a)


