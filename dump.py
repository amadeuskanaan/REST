
'#####03 '
        # FSL
        # if not os.path.isfile('REST_DISCO_MOCO1.nii.gz'):
        #     os.system('fslsplit ../REST_DISCO.nii.gz -t')
        #     for i in xrange(0,417):
        #        frame = '{0:0>4}'.format(i)
        #        os.system('flirt -in vol%s -ref REST_DICSO_mean.nii.gz -out flirt_1_%s.nii.gz'%(frame,frame))
        #     os.system('fslmerge -t REST_DISCO_MOCO1.nii.gz flirt_1*')
        #     os.system('fslmaths REST_DISCO_MOCO1 -Tmean REST_DISCO_MOCO1_mean.nii.gz')
        #     os.system('rm -rf flirt_1*')


        # #2 MOCO step 2
        # if not os.path.isfile('REST_DISCO_MOCO2.nii.gz'):
        #     for i in xrange(0,417):
        #        frame = '{0:0>4}'.format(i)
        #        os.system('flirt -in vol%s -ref REST_DISCO_MOCO1_mean.nii.gz -out flirt_2_%s.nii.gz -omat omat_2_%s.mat'%(frame, frame, frame))
        #     os.system('fslmerge -t REST_DISCO_MOCO2.nii.gz flirt_2*')
        #     os.system('fslmaths REST_DISCO_MOCO2.nii.gz -Tmean REST_DICSO_MOCO2_mean.nii.gz')
        #     os.system('cp REST_DISCO_MOCO2.nii.gz ../REST_DISCO_MOCO2.nii.gz')
        #     os.system('rm -rf flirt_2* vol*')


'#####04 '

            #    os.system('applywarp '
            #               '-i %s/CONVERT_XFM/vol%s '   ##### vols are splits of  os.path.join(subject_dir, 'FUNC_PPROC/REST_DROP_RPI.nii.gz')
            #               '-o applywarp%s.nii.gz '
            #               '-r FUNC2ANAT_STEP_2_BBR.nii.gz '
            #               '--postmat=MATS_UNIFIED/UNIMAT_%s.mat '
            #               '--rel '
            #               '-w %s/CONVERT_XFM/warp%s.nii.gz'
            #               %(pproc_dir, frame, frame, frame, pproc_dir, frame))

        #     os.system('fslmerge -t REST_PPROC_ANAT2mm.nii.gz applywarp*')
        #     os.system('fslmerge -t REST_UNIFIED_ANAT2mm.nii.gz applywarp*')