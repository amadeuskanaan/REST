__author__ = 'kanaan 11.12.2015'

import os
from variables.subject_list import *
from utilities.nuisance import *
from utilities.qc  import *
import pandas as pd
from nipype.algorithms.misc import TSNR
from reportlab.pdfgen import canvas
from reportlab.lib.units import cm, mm, inch, pica
from nilearn import plotting

def get_distributions(population, workspace_dir):
    fd_means           = []
    tsnr_brain_medians = []
    tsnr_medians_first = []

    print 'Grabbing FD and TSNR population distribution'
    for subject in population:
        # print subject
        subject_dir       = os.path.join(workspace_dir, subject)
        qcdir             = os.path.join(subject_dir,'QUALITY_CONTROL')
        mkdir_path(qcdir)
        os.chdir(qcdir)
        print subject

        #fd
        movpar            = os.path.join(subject_dir,'FUNC_PPROC/MOTION_CORRECTION/REST_DISCO_MOCO_2.1D')
        fd = calc_FD_power(movpar)
        fd_means.append(np.mean(np.loadtxt(fd)))

        #tsnr
        func_native       = os.path.join(subject_dir, 'FUNC_PPROC/REST_PPROC_NATIVE_BRAIN.nii.gz')
        func_native_mask  = os.path.join(subject_dir, 'FUNC_PPROC/REST_PPROC_NATIVE_BRAIN_mask_ero.nii.gz')
        func_native_first =  os.path.join(subject_dir, 'FUNC_TRANSFORM/NATIVE_FUNC_FIRST.nii.gz')

        # if not os.path.join(qcdir, 'REST_PPROC_NATIVE_BRAIN_tsnr.nii.gz'):
        print 'tsnr'
        if not os.path.isfile(os.path.join(qcdir, 'REST_PPROC_NATIVE_BRAIN_tsnr.nii.gz')):
            tsnr = TSNR()
            tsnr.inputs.in_file = func_native
            tsnr.run()

        tsnr_data    = nb.load('REST_PPROC_NATIVE_BRAIN_tsnr.nii.gz').get_data()
        nan_mask     = np.logical_not(np.isnan(tsnr_data))
        brain_mask   = nb.load(func_native_mask).get_data() > 0
        first_mask   = nb.load(func_native_first).get_data() > 0

        tsnr_brain_median = np.median(tsnr_data[np.logical_and(nan_mask, brain_mask)])
        tsnr_first_median = np.median(tsnr_data[np.logical_and(nan_mask, first_mask)])

        print tsnr_brain_median, tsnr_first_median

        tsnr_brain_medians.append(tsnr_brain_median)
        tsnr_medians_first.append(tsnr_first_median)


    np.savetxt(os.path.join(workspace_dir, 'population_fd_distributions.txt'), fd_means)
    np.savetxt(os.path.join(workspace_dir, 'population_tsnr_distributions.txt'), tsnr_brain_medians)
    np.savetxt(os.path.join(workspace_dir, 'population_tsnr_first_distributions.txt'), tsnr_medians_first)

def make_quality_control_reports(population, workspace_dir):

    tsnr_dist = []
    fd_dist   = []
    for subject in population:

        print '###############################################################################'
        print ' Running Quality Control for subject %s' %subject
        print ''

        #input
        subject_dir       = os.path.join(workspace_dir, subject)

        func_native       = os.path.join(subject_dir, 'FUNC_PPROC/REST_PPROC_NATIVE_BRAIN.nii.gz')
        func_native_mask  = os.path.join(subject_dir, 'FUNC_PPROC/REST_PPROC_NATIVE_BRAIN_mask_ero.nii.gz')
        func_native_mean  = os.path.join(subject_dir, 'FUNC_PPROC/REST_PPROC_NATIVE_BRAIN_mean.nii.gz')
        func_native_gm    = os.path.join(subject_dir, 'FUNC_TRANSFORM/NATIVE_FUNC_GM.nii.gz')
        func_native_wm    = os.path.join(subject_dir, 'FUNC_TRANSFORM/NATIVE_FUNC_WM.nii.gz')
        func_native_csf   = os.path.join(subject_dir, 'FUNC_TRANSFORM/NATIVE_FUNC_CSF.nii.gz')
        func_native_first =  os.path.join(subject_dir, 'FUNC_TRANSFORM/NATIVE_FUNC_FIRST.nii.gz')

        func_mni          = os.path.join(subject_dir, 'FUNC_TRANSFORM/REST_PPROC_MNI2mm_BRAIN.nii.gz')
        func_mni_mask     = os.path.join(subject_dir, 'FUNC_TRANSFORM/REST_PPROC_MNI2mm_BRAIN_mask_ero.nii.gz')
        func_mni_mean     = os.path.join(subject_dir, 'FUNC_TRANSFORM/REST_PPROC_MNI2mm_BRAIN_mean.nii.gz')
        func_mni_gm       = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_GM.nii.gz')
        func_mni_wm       = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_WM.nii.gz')
        func_mni_csf      = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_CSF.nii.gz')

        residual_compor   = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_compcor_friston_bp_fwhm.nii.gz')
        residual_wmcsf    = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_wmcsf_friston_bp_fwhm.nii.gz')
        residual_global   = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_detrend_global_wmcsf_friston_bp_fwhm.nii.gz')

        aroma_compcor     = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_compcor_friston_bp.nii.gz')
        aroma_wmcsf       = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_wmcsf_friston_bp.nii.gz')
        aroma_global      = os.path.join(subject_dir, 'FUNC_DENOISE/RESIDUAL_MNI2mm_FWHM_AROMA_detrend_global_wmcsf_friston_bp.nii.gz')

        movpar            = os.path.join(subject_dir,'FUNC_PPROC/MOTION_CORRECTION/REST_DISCO_MOCO_2.1D')

        qcdir =  os.path.join(subject_dir,'QUALITY_CONTROL')
        mkdir_path(qcdir)
        os.chdir(qcdir)

        print 'Creating Quality Control Report'

        ################################################################################################################
        print '1. Calculating FD'
        fd = os.path.join(qcdir, 'FD.1D')
        population_FDs = np.genfromtxt(os.path.join(workspace_dir, 'population_fd_distributions.txt'))
        if not os.path.isfile('plot_fd_qc.png'):
           plot_FD(fd, population_FDs, subject,figsize = (8.3,8.3))

        ################################################################################################################
        print '2. Grabbing accepted frames (FD<0.2mm) '
        fd1d = np.loadtxt(fd)
        in_frames = []
        for frame, fd in enumerate(fd1d):
            if fd < 0.2:
                in_frames.append(frame)
        print '    ...Subject has %s of 417 good frames'%len(in_frames)
        np.save('in_frames', str(in_frames).replace(" ",""))
        if len(in_frames) > 150:
            in_frames_100 = in_frames[0:130]
            np.save('in_frames_100', str(in_frames_100).replace(" ",""))
        fd_dist.append(np.mean(fd))

        ################################################################################################################
        print '3. Calculating DVARS'
        dvars = os.path.join(qcdir, 'DVARS.npy')
        if not os.path.isfile(dvars):
            calc_DVARS(func_native, func_native_mask)

        ###############################################################################################################
        print '4. TSNR'
        #get TSNR median for whole brain
        tsnr_data   = nb.load('./REST_PPROC_NATIVE_BRAIN_tsnr.nii.gz').get_data()
        nan_mask    = np.logical_not(np.isnan(tsnr_data))
        mask        = nb.load(func_native_mask).get_data() > 0
        tsnr_median = np.median(tsnr_data[np.logical_and(nan_mask, mask)])
        tsnr_dist.append(tsnr_median)

        if not os.path.isfile(' plot_tsnr_mosaic.png'):
            plot_mosaic(nifti_file  = './REST_PPROC_NATIVE_BRAIN_tsnr.nii.gz',
                        output_name = 'plot_tsnr_mosaic.pdf',
                        title="tSNR",
                        overlay_mask = None,
                        figsize=(8.3, 11.7))

            os.system('convert -density 300 -trim plot_tsnr_mosaic.pdf -quality 300 -sharpen 0x1.0 -crop  2506x2570+1+470 plot_tsnr_mosaic.png')

        # whole brain plots
        if not os.path.isfile('plot_tsnr_brain_distribution.png'):
            plot_distrbution_of_values(main_file =  './REST_PPROC_NATIVE_BRAIN_tsnr.nii.gz',
                                       mask_file = func_native_mask,
                                       xlabel = "%s Whole Brain tSNR distribution" % subject,
                                       distribution=  np.genfromtxt(os.path.join(workspace_dir, 'population_tsnr_distributions.txt')),
                                       xlabel2= "%s median whole brain tSNR with respect to all subjects"%subject,
                                       figsize=(11.7,8.3),
                                       outname = 'plot_tsnr_brain_distribution.png')

        d = plotting.plot_epi(func_native_mean, cmap='gray', display_mode='y', black_bg= 1, annotate=0)
        d.add_contours(func_native_mask, colors='r')
        d.savefig('plot_func_native_mask_y.png', dpi = 130)

        d = plotting.plot_epi(func_native_mean, cmap='gray', display_mode='x', black_bg= 1, annotate=0)
        d.add_contours(func_native_mask, colors='r')
        d.savefig('plot_func_native_mask_x.png', dpi = 110)

        d = plotting.plot_epi(func_native_mean, cmap='gray', display_mode='z', black_bg= 1, annotate=0)
        d.add_contours(func_native_mask, colors='r')
        d.savefig('plot_func_native_mask_z.png', dpi = 130)

        # subcortical plots
        if not os.path.isfile('plot_tsnr_first_distribution.png'):
            plot_distrbution_of_values(main_file    =  './REST_PPROC_NATIVE_BRAIN_tsnr.nii.gz',
                                       mask_file    = func_native_first,
                                       xlabel       = "%s Subcortical tSNR distribution" % subject,
                                       distribution =  np.genfromtxt(os.path.join(workspace_dir, 'population_tsnr_first_distributions.txt')),
                                       xlabel2      = "%s median subcortical tSNR with respect to all subjects"%subject,
                                       figsize      =(11.7,8.3),
                                       outname = 'plot_tsnr_first_distribution.png')

            plot_3d_overlay(underlay_file=func_native_mean, overlay_file=func_native_first, out_filename='plot_subcortical_mask.png', dpi = 215)

        ################################################################################################################
        print '5. Plot FUNC 2 MNI Registration'

        d = plotting.plot_epi( mni_brain_2mm, cmap='gray', display_mode='x',  annotate=0, black_bg=0)
        d.add_contours(func_mni_mean, colors='r')
        d.savefig('plot_mnireg_x.png', dpi = 110)

        d = plotting.plot_epi( mni_brain_2mm, cmap='gray', display_mode='y',  annotate=0, black_bg=0)
        d.add_contours(func_mni_mean, colors='r')
        d.savefig('plot_mnireg_y.png', dpi = 130)

        d = plotting.plot_epi( mni_brain_2mm, cmap='gray', display_mode='z',  annotate=0, black_bg=1)
        d.add_contours(func_mni_mean, colors='r')
        d.savefig('plot_mnireg_z.png', dpi = 130)


        ################################################################################################################
        print '6. Plot Func GM'
        if not os.path.isfile('plot_func_native_gm.png'):
            plot_mosaic(nifti_file  = func_native_mean,
                        output_name = 'plot_func_native_gm.pdf',
                        title="Func Gray Matter",
                        overlay_mask = func_native_gm,
                        figsize=(8.3, 11.7))
            os.system('convert -density 300 -trim plot_func_native_gm.pdf -quality 300 -sharpen 0x1.0 -crop 2506x2570+1+470 plot_func_native_gm.png')

        #
        ################################################################################################################
        print '7. Plot Nuisance Residuals'

        mni_wm_sig  = os.path.join(subject_dir, 'FUNC_DENOISE/TISSUE_SIGNALS_MNI/NUISANCE_SIGNALS_WM.npy')
        mni_gm_sig  = os.path.join(subject_dir, 'FUNC_DENOISE/TISSUE_SIGNALS_MNI/NUISANCE_SIGNALS_GM.npy')
        mni_csf_sig = os.path.join(subject_dir, 'FUNC_DENOISE/TISSUE_SIGNALS_MNI/NUISANCE_SIGNALS_CSF.npy')

        if not os.path.isfile('func_mni_detrend.nii.gz'):
           calc_residuals(subject = func_mni, selector = nuisance_detrend,wm_sig_file = mni_wm_sig,
                          csf_sig_file= mni_csf_sig, gm_sig_file = mni_gm_sig,motion_file = movpar,
                          compcor_ncomponents = 0)
           os.system('mv residual.nii.gz func_mni_detrend.nii.gz')
           os.system('rm -rf  quadratic_constant_linear.csv nuisance_regressors.mat')

        if not os.path.isfile('plot_nuisance.png'):
            plot_nuisance_residuals(mov_params                       = movpar,
                                    fd1d                             = fd1d,
                                    func_preprocessed                = func_native,
                                    func_preprocessed_mask           = func_native_mask,
                                    dvars                            = dvars,
                                    func_gm                          = func_mni_gm,
                                    residuals_dt                     = 'func_mni_detrend.nii.gz',
                                    residuals_cc                     = residual_compor ,
                                    residuals_gl                     = residual_global,
                                    aroma_cc                         = aroma_compcor ,
                                    aroma_gl                         = aroma_global,
                                    out_name                         = 'plot_nuisance.png',
                                    figsize = (8.3,8.3))

        ################################################################################################################
        print '8. Calculating Motion statistics and saving as csv'
        df = pd.DataFrame(index = [subject], columns = ['FD', 'FD_in','FD_topQua', 'FD_max', 'FD_RMS', 'DVARS', 'TSNR'])
        df.loc[subject]['FD']        = str(np.round(np.mean(fd1d), 3))
        df.loc[subject]['FD_in']     = str(np.round(len(in_frames), 3))
        quat = int(len(fd1d)/4)
        df.loc[subject]['FD_topQua'] = str(np.round(np.mean(np.sort(fd1d)[::-1][:quat]), 3))
        df.loc[subject]['FD_max']    = str(np.round(np.max(fd1d), 3))
        df.loc[subject]['FD_RMS']    = str(np.round(np.sqrt(np.mean(fd1d)), 3))
        df.loc[subject]['DVARS']     = str(np.round(np.mean(np.load(dvars)), 3))
        df.loc[subject]['TSNR']      = str(np.round(tsnr_median, 3))
        df.to_csv ('quality_paramters.csv')

        ################################################################################################################

        print 'Creating QC REPORT'

        report = canvas.Canvas('QC_REPORT.pdf', pagesize=(1280 *1.9, 1556*1.9))
        report.setFont("Helvetica", 100)
        report.drawString(inch*15, inch*39.3, '%s'%subject)
        report.setFont("Helvetica", 70)
        report.drawString(inch*15, inch*1, 'tSNR')
        report.drawImage('plot_tsnr_mosaic.png', inch*1.5, inch*2.6)
        report.showPage()

        report.setFont("Helvetica", 30)
        report.drawString(inch*31, inch*40, '%s'%subject)
        report.setFont("Helvetica", 70)
        report.drawString(inch*12, inch*1, 'Func Native Gray Matter')
        report.drawImage('plot_func_native_gm.png', inch*1.5, inch*2.6)
        report.showPage()

        report.setFont("Helvetica", 30)
        report.drawString(inch*31, inch*40, '%s'%subject)
        report.setFont("Helvetica", 50)
        report.drawString(inch*10, inch*38.5, 'Func to MNI2mm Nonlinear Registration')
        report.drawImage('plot_mnireg_x.png', inch*2.5, inch*34.5)
        report.drawImage('plot_mnireg_y.png', inch*2.5, inch*30.5)

        report.setFont("Helvetica", 30)
        report.drawString(inch*31, inch*40, '%s'%subject)
        report.setFont("Helvetica", 50)
        report.drawString(inch*13, inch*29, 'Func Native Brain Mask')
        report.drawImage('plot_func_native_mask_y.png', inch*2.5, inch*21)
        report.drawImage('plot_func_native_mask_x.png', inch*2.5, inch*25)
        report.drawImage('plot_tsnr_brain_distribution.png', inch*3.5, inch*1.7)
        report.showPage()

        report.setFont("Helvetica", 30)
        report.drawString(inch*31, inch*40, '%s'%subject)
        report.setFont("Helvetica", 50)
        report.drawString(inch*12, inch*38.5, 'Func FIRST Subcortical Mask')
        report.drawImage('plot_subcortical_mask.png', inch*7, inch*20)
        report.drawImage('plot_tsnr_first_distribution.png', inch*3, inch*2)
        report.showPage()

        report.setFont("Helvetica", 30)
        report.drawString(inch*31, inch*40, '%s'%subject)
        report.drawImage('plot_nuisance.png', inch*3, inch*2)
        report.drawImage('plot_fd_qc.png', inch*3, inch*20)
        report.showPage()
        report.save()

def create_population_qc_spearsheet(population, workspace):
    print ''

# get_distributions(controls_a + patients_a , workspace_a)
# get_distributions(patients_b , workspace_b)
make_quality_control_reports(['BM8X'], workspace_a)