import os


#def plot_surf(working_dir):
from surfer import Brain, io
#ecm = os.path.join(working_dir, 'STATISTICS/ECM/%s'%string)
ecm = '/scr/sambesi4/workspace/project_REST/study_a/STATISTICS/SCA/FIRST/RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm_scrubbed/randomise_baseline_tfce_corrp_tstat2.nii.gz'
reg_file = '/afs/cbs.mpg.de/software/freesurfer/5.3.0/ubuntu-precise-amd64/average/mni152.register.dat'
surf_data_lh = io.project_volume_data(ecm, "lh", reg_file)
brain = Brain("fsaverage", "lh", "pial", views=['lat', 'med'], background="black")
brain.add_data(surf_data_lh, min = 0, max = 1,   colormap="jet", hemi='lh')
# brain.save_image('/SCR/ECM_imgs/all_mean.png')
