__author__ = 'kanaan 06.01.2016'

import os
from utilities.utils import *
from variables.subject_list import *
from nilearn import input_data
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

columns =['cau_put', 'cau_pal', 'cau_tha', 'cau_stn', 'cau_sn', 'cau_acc', 'cau_ins',
                     'put_pal', 'put_tha', 'put_stn', 'put_sn', 'put_acc', 'put_ins',
                                'pal_tha', 'pal_stn', 'pal_sn', 'pal_acc', 'pal_ins',
         ]


def run_rsfc(population, workspace_dir, population_name):

    df = pd.DataFrame(index = [population], columns= columns)

    for resid in resid_strings:
        for subject in population:
            print '#############################################'
            print 'Running RSFC for Subject:%s'%subject
            print ''

            # input/ouput
            subject_dir    = os.path.join(workspace_dir, subject)
            rsfc_dir       = os.path.join(subject_dir, 'RSFC/%s'%resid)
            mkdir_path(rsfc_dir)
            os.chdir(rsfc_dir)

            # grab func and mask images
            if not resid[-8:] == 'scrubbed':
                pproc    = os.path.join(subject_dir, 'FUNC_DENOISE/%s.nii.gz'%resid)
            elif resid[-8:] == 'scrubbed':
                pproc    = os.path.join(subject_dir, 'FUNC_DENOISE_SCRUB/%s.nii.gz'%resid)

            caudate  = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_CAUDATE_L.nii.gz')
            putamen  = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_PUTAMEN_L.nii.gz')
            pallidum = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_PALLIDUM_L.nii.gz')
            thalamus = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_THALAMUS_L.nii.gz')
            insula   = os.path.join(subject_dir, 'FUNC_TRANSFORM/MNI2mm_FUNC_INSULA_LA.nii.gz')
            acc      = '/SCR/ROI/ACCsphere_5mm.nii.gz'
            atag_stn = '/scr/sambesi1/workspace/Projects/REST/atlases/ATAG_2mm_STN_bin.nii.gz'
            atag_sn  = '/scr/sambesi1/workspace/Projects/REST/atlases/ATAG_2mm_SN_bin.nii.gz'


            print 'Grab timeseries and calculating correlations for', resid
            if not os.path.isfile('caudate.txt'):
                print '.........extracting caudate'
                np.savetxt('caudate.txt',  input_data.NiftiLabelsMasker(labels_img= caudate,  standardize=True).fit_transform(pproc))
            if not os.path.isfile('putamen.txt'):
                print '.........extracting putamen'
                np.savetxt('putamen.txt',  input_data.NiftiLabelsMasker(labels_img= putamen,  standardize=True).fit_transform(pproc))
            if not os.path.isfile('pallidum.txt'):
                print '.........extracting pallidum'
                np.savetxt('pallidum.txt', input_data.NiftiLabelsMasker(labels_img= pallidum, standardize=True).fit_transform(pproc))
            if not os.path.isfile('thalamus.txt'):
                print '.........extracting thalamus'
                np.savetxt('thalamus.txt', input_data.NiftiLabelsMasker(labels_img= thalamus, standardize=True).fit_transform(pproc))
            if not os.path.isfile('stn.txt'):
                print '.........extracting stn'
                np.savetxt('stn.txt', input_data.NiftiLabelsMasker(labels_img= atag_stn, standardize=True).fit_transform(pproc))
            if not os.path.isfile('sn.txt'):
                print '.........extracting sn'
                np.savetxt('sn.txt', input_data.NiftiLabelsMasker(labels_img= atag_sn, standardize=True).fit_transform(pproc))
            if not os.path.isfile('acc.txt'):
                print '.........extracting acc'
                np.savetxt('acc.txt', input_data.NiftiLabelsMasker(labels_img= acc, standardize=True).fit_transform(pproc))
            if not os.path.isfile('insula.txt'):
                print '.........extracting insula'
                np.savetxt('insula.txt', input_data.NiftiLabelsMasker(labels_img= insula, standardize=True).fit_transform(pproc))

            timeseries_cau = np.genfromtxt('caudate.txt')
            timeseries_put = np.genfromtxt('putamen.txt')
            timeseries_pal = np.genfromtxt('pallidum.txt')
            timeseries_tha = np.genfromtxt('thalamus.txt')
            timeseries_stn = np.genfromtxt('stn.txt')
            timeseries_sn  = np.genfromtxt('sn.txt')
            timeseries_acc = np.genfromtxt('acc.txt')
            timeseries_ins = np.genfromtxt('insula.txt')

            # CAUDATE
            df.loc[subject]['cau_put']  = float(pearsonr(timeseries_cau, timeseries_put)[0])
            df.loc[subject]['cau_pal']  = float(pearsonr(timeseries_cau, timeseries_pal)[0])
            df.loc[subject]['cau_tha']  = float(pearsonr(timeseries_cau, timeseries_tha)[0])
            df.loc[subject]['cau_stn']  = float(pearsonr(timeseries_cau, timeseries_stn)[0])
            df.loc[subject]['cau_sn']   = float(pearsonr(timeseries_cau, timeseries_sn)[0])
            df.loc[subject]['cau_acc']  = float(pearsonr(timeseries_cau, timeseries_acc)[0])
            df.loc[subject]['cau_ins']  = float(pearsonr(timeseries_cau, timeseries_ins)[0])

            # PUTAMEN
            df.loc[subject]['put_pal']  = float(pearsonr(timeseries_put, timeseries_pal)[0])
            df.loc[subject]['put_tha']  = float(pearsonr(timeseries_put, timeseries_tha)[0])
            df.loc[subject]['put_stn']  = float(pearsonr(timeseries_put, timeseries_stn)[0])
            df.loc[subject]['put_sn']   = float(pearsonr(timeseries_put, timeseries_sn)[0])
            df.loc[subject]['put_acc']  = float(pearsonr(timeseries_put, timeseries_acc)[0])
            df.loc[subject]['put_ins']  = float(pearsonr(timeseries_put, timeseries_ins)[0])

            # PALLIDUM
            df.loc[subject]['pal_tha']  = float(pearsonr(timeseries_pal, timeseries_tha)[0])
            df.loc[subject]['pal_stn']  = float(pearsonr(timeseries_pal, timeseries_stn)[0])
            df.loc[subject]['pal_sn']   = float(pearsonr(timeseries_pal, timeseries_sn)[0])
            df.loc[subject]['pal_acc']  = float(pearsonr(timeseries_pal, timeseries_acc)[0])
            df.loc[subject]['pal_ins']  = float(pearsonr(timeseries_pal, timeseries_ins)[0])

        mkdir_path(os.path.join(workspace_dir, 'STATISTICS/RSFC/%s'%population_name))
        df.to_csv(os.path.join(workspace_dir, 'STATISTICS/RSFC/%s/RSFC_%s.csv'%(population_name, resid)))


# run_rsfc(controls_a, workspace_a, 'controls')
run_rsfc(patients_a, workspace_a, 'patients')
# run_rsfc(patients_b, workspace_b, 'patients')








