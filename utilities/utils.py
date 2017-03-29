__author__ = 'kanaan'
from nipype.pipeline.engine import Workflow, Node
import nipype.interfaces.utility as util
from nipype.interfaces.afni import preprocess
import nipype.interfaces.freesurfer as fs

import string
valid_chars = '-_.() %s%s' %(string.ascii_letters, string.digits)


def grabber_util(name):
    '''
    Method to grab data from from out_dir
    inputnode
        input.base_dir
        input.subject_id
        input.folder_name
        input.file_name
        input.pipeline_name
    outputnode
        output.file
    '''

    def grab_filepath(base_dir,subject_id, folder_name, file_name_string):
        import os
        for file in os.listdir(os.path.join(base_dir,
                                             subject_id,
                                             folder_name)):
            if file_name_string in file:
                file_path =os.path.join(base_dir,
                                        subject_id,
                                        folder_name,
                                        file)
                return file_path


    flow      = Workflow(name=name)
    inputnode = Node(util.IdentityInterface(fields = ['subject_id','folder_name', 'file_name_string', 'base_dir']),
                        name = 'inputnode')
    outputnode = Node(util.IdentityInterface(fields=['out_file']),
                        name = 'outputnode')


    grabber  = Node(util.Function(input_names      = ['base_dir', 'subject_id', 'folder_name', 'file_name_string'],
                                  output_names     = ['out_file'],
                                  function         = grab_filepath),
                                  name             = 'grabber')

    flow.connect(inputnode   , 'subject_id'       ,   grabber,     'subject_id'         )
    flow.connect(inputnode   , 'base_dir'         ,   grabber,     'base_dir'           )
    flow.connect(inputnode   , 'folder_name'      ,   grabber,     'folder_name'        )
    flow.connect(inputnode   , 'file_name_string' ,   grabber,     'file_name_string'   )
    flow.connect(grabber     , 'out_file'         ,   outputnode,  'out_file'           )


    return flow

def return_list(file_1,file_2):
    list = [file_1,file_2]
    return list

def locate(string, directory):
        import os
        for file in os.listdir(directory):
            x=[]
            if string in file:
                x = os.path.join(directory,file)
                return x


import string
valid_chars = '-_.() %s%s' %(string.ascii_letters, string.digits)


def mkdir_path(path):
    import os
    import errno
    import string
    import subprocess

    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise






def calc_friston_twenty_four(mov_par):
    import numpy as np
    import os
    twenty_four   = None

    six           = np.genfromtxt(mov_par)
    six_squared   = six**2

    twenty_four   = np.concatenate((six,six_squared), axis=1)

    six_roll      = np.roll(six, 1, axis=0)
    six_roll[0]   = 0

    twenty_four   = np.concatenate((twenty_four, six_roll), axis=1)

    six_roll_squ  = six_roll**2

    twenty_four   = np.concatenate((twenty_four, six_roll_squ), axis=1)
    updated_mov   = os.path.join(os.getcwd(), 'FRISTON_24.1D')
    np.savetxt(updated_mov, twenty_four, fmt='%0.8f', delimiter=' ')

    return updated_mov


def safe_shape(*vol_data):
    """
    Checks if the volume (first three dimensions) of multiple ndarrays
    are the same shape.

    Parameters
    ----------
    vol_data0, vol_data1, ..., vol_datan : ndarray
        Volumes to check

    Returns
    -------
    same_volume : bool
        True only if all volumes have the same shape.
    """
    same_volume = True

    first_vol_shape = vol_data[0].shape[:3]
    for vol in vol_data[1:]:
        same_volume &= (first_vol_shape == vol.shape[:3])

    return same_volume