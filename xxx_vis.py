import os


#def plot_surf(working_dir):
from surfer import Brain, io
#ecm = os.path.join(working_dir, 'STATISTICS/ECM/%s'%string)
ecm = '/scr/sambesi4/workspace/project_REST/study_a/STATISTICS/SCA/FIRST/RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm_scrubbed/randomise_baseline_tfce_corrp_tstat2.nii.gz'
ecm = '/scr/sambesi3/workspace/project_iron/study_a/STATISTICS/randomise_qsm_tfce_corrp_tstat1.nii.gz'
reg_file = '/afs/cbs.mpg.de/software/freesurfer/5.3.0/ubuntu-precise-amd64/average/mni152.register.dat'
surf_data_lh = io.project_volume_data(ecm, "lh", reg_file)
brain = Brain("fsaverage", "lh", "pial", views=['lat', 'med'], background="black")
brain.add_data(surf_data_lh, min = 0, max = 1,   colormap="jet", hemi='lh')
# brain.save_image('/SCR/ECM_imgs/all_mean.png')


from surfer import Brain, io

x = '/scr/sambesi3/workspace/project_iron/study_a/HCTT/SEGMENTATION/FREESURFER/ins.nii.gz'
reg_file = '/afs/cbs.mpg.de/software/freesurfer/5.3.0/ubuntu-precise-amd64/average/mni152.register.dat'
surf_data_lh = io.project_volume_data(x, "lh", reg_file)
surf_data_lh = io.project_volume_data(x, "rh", reg_file)
brain = Brain("fsaverage", "split", "pial", views=['lat', 'med'], background="white")
brain.add_data(surf_data_lh, min = 0, max = 1, thresh = 0.9,   colormap="jet", hemi='lh')
brain.add_data(surf_data_lh, min = 0, max = 1, thresh = 0.9,   colormap="jet", hemi='rh')