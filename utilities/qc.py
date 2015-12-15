# see https://github.com/FCP-INDI/C-PAC
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
import os
import math
import time
import seaborn as sns
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import FigureCanvasPdf as FigureCanvas


def calc_FD_power(motion_pars):
    '''
    Method to calculate FD based on (Power, 2012)
    '''
    import os
    import numpy as np

    fd_out       =  os.path.join(os.getcwd(), 'FD.1D')
    lines        =  open(motion_pars, 'r').readlines()
    rows         = [[float(x) for x in line.split()] for line in lines]
    cols         = np.array([list(col) for col in zip(*rows)])
    translations = np.transpose(np.abs(np.diff(cols[0:3, :])))
    rotations    = np.transpose(np.abs(np.diff(cols[3:6, :])))
    FD_power     = np.sum(translations, axis = 1) + (50*3.141/180)*np.sum(rotations, axis =1)
    #FD is zero for the first time point
    FD_power = np.insert(FD_power, 0, 0)

    np.savetxt(fd_out, FD_power)

    return fd_out

def calc_DVARS(rest, mask):
    '''
    Method to calculate DVARS according to (Power, 2012)
    CPAC-0.3.8 implenentation
    '''
    import numpy as np
    import nibabel as nb
    import os

    dvars_out    = os.path.join(os.getcwd(), 'DVARS.npy')
    rest_data    = nb.load(rest).get_data().astype(np.float32)
    mask_data    = nb.load(mask).get_data().astype(np.bool)
    #square of relative intensity value for each voxel across every timepoint
    data         = np.square(np.diff(rest_data, axis = 3))
    #applying mask, getting the data in the brain only
    data         = data[mask_data]
    #square root and mean across all timepoints inside mask
    DVARS        = np.sqrt(np.mean(data, axis=0))

    np.save(dvars_out, DVARS)
    return dvars_out

def gen_realignment_params(realignment_parameters_file):
    data = np.loadtxt(realignment_parameters_file)
    data_t = data.T
    x = data_t[0]
    y = data_t[1]
    z = data_t[2]
    for i in range(3, 6):
        for j in range(len(data_t[i])):
            data_t[i][j] = math.degrees(data_t[i][j])
    roll = data_t[3]
    pitch= data_t[4]
    yaw = data_t[5]
    return x,y,z,roll, pitch, yaw

def timeseries(rest, grey):
    import numpy as np
    import nibabel as nib
    import os
    rest_data = nib.load(rest).get_data().astype(np.float32)
    gm_mask = nib.load(grey).get_data().astype('bool')
    rest_gm = rest_data[gm_mask]

    return rest_gm


def plot_vline(cur_val, label, ax):
    ax.axvline(cur_val)
    ylim = ax.get_ylim()
    vloc = (ylim[0] + ylim[1]) / 2.0
    xlim = ax.get_xlim()
    pad = (xlim[0] + xlim[1]) / 100.0
    ax.text(cur_val - pad, vloc, label, color="blue", rotation=90, verticalalignment='center', horizontalalignment='right')


def plot_FD(fd1d, mean_FD_distribution, subject,figsize = (8.3,8.3)):

    threshold = 0.2
    FD_power      = np.genfromtxt(fd1d)
    meanFD        = np.round(np.mean(FD_power), decimals = 2)
    rmsFD         = np.sqrt(np.mean(FD_power))
    count         = np.int(FD_power[FD_power>threshold].size)
    percentFD     = (count*100/(len(FD_power)+1))

    fig = plt.figure(figsize=figsize)

    fig.subplots_adjust(wspace=0.3)
    fig.set_size_inches(12, 8)

    grid = GridSpec(2, 4)

    ax = plt.subplot(grid[0,:-1])
    ax.plot(FD_power)
    ylim = ax.get_ylim()
    ax.set_xlim((0, len(FD_power)))
    ax.set_ylabel("%s Frame Displacement [mm]" %subject)
    ax.set_xlabel("Frame #")
    ax.text(20, (ylim[1] - ylim[1]*0.05), 'FD mean = %s'%meanFD, va='center', size = 18, color = 'r')
    ax.text(20, (ylim[1] - ylim[1]*0.2), '%s Frames (%s%%) > threshold  '%(count, percentFD), va='center', size = 18, color = 'r')
    #ax.text(20, (ylim[1] - ylim[1]*0.15), 'are above the threshold '%, va='center', size = 18, color = 'r')
    plt.axhline(threshold, linestyle='dashed', linewidth=2)#,color='r')

    ax = plt.subplot(grid[0,-1])
    sns.distplot(FD_power, vertical = True, ax = ax)
    ax.set_ylim(ylim)

    ax= plt.subplot(grid[1,:])
    sns.distplot(mean_FD_distribution, ax=ax)
    ax.set_xlabel("%s Mean Frame Dispalcement (over all subjects) [mm]"%subject)
    label = "MeanFD = %g"%meanFD
    plot_vline(meanFD, label, ax=ax)

    png_name = str(os.path.join(os.getcwd(), 'plot_fd_qc.png'))
    plt.savefig(png_name, dpi=190, bbox_inches='tight')



def plot_nuisance_residuals(mov_params,
                            fd1d,
                            func_preprocessed,
                            func_preprocessed_mask,
                            dvars,
                            func_gm,
                            residuals_dt,
                            residuals_cc,
                            residuals_gl,
                            aroma_cc,
                            aroma_gl,
                            out_name,
                            figsize = (8.3,8.3)
                            ):

    sns.set_style('darkgrid')

    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.05)
    fig.set_size_inches(12, 8)

    # open figure
    fig = plt.figure()
    fig.subplots_adjust(wspace=0.05)
    fig.set_size_inches(12, 8)

    # 1 plot translations
    movement = gen_realignment_params(mov_params)
    ax1 = plt.subplot2grid((8,2), (0,0),  colspan = 1, rowspan =1)
    ax1.yaxis.set_tick_params(labelsize=5)
    ax1.plot(movement[0])
    ax1.plot(movement[1])
    ax1.plot(movement[2])
    ax1.set_xticklabels([])
    ax1.grid(True)
    ax1.set_color_cycle(['red', 'green', 'blue'])
    ax1.set_ylabel('trans')

    # 2 plot rotations
    ax2 = plt.subplot2grid((8,2), (0,1),  colspan = 1, rowspan =1)
    ax2.set_color_cycle(['red', 'green', 'blue'])
    ax2.plot(movement[3])
    ax2.plot(movement[4])
    ax2.plot(movement[5])
    ax2.yaxis.set_tick_params(labelsize=5)
    ax2.set_xticklabels([])
    ax2.grid(True)
    ax2.set_ylabel('rot')

    # 3 plot FD
    FD_power = np.loadtxt(fd1d)
    ax3 = plt.subplot2grid((8,2), (1, 0),  colspan = 2, rowspan =1)
    #ax3.axes.get_yaxis().set_visible(True)
    ax3.yaxis.set_tick_params(labelsize=5)
    ax3.set_xticklabels([])
    ax3.grid(True)
    ax3.plot(FD_power)
    ax3.set_xlim([0, 422])
    plt.axhline(0.2, linestyle='dashed', linewidth=2, color='r')
    ax3.set_ylabel('FD')
    #ax3.yaxis.set_major_locator(FixedLocator((ax3.get_ylim())))

    # 4 plot DVARS
    # DVARS = calc_DVARS(func_preprocessed, func_preprocessed_mask)
    # DVARS = calc_dvars(func_preprocessed, output_all=False, interp="fraction")
    dv = np.load(dvars)
    ax4 = plt.subplot2grid((8,2), (2, 0),  colspan = 2, rowspan =1)
    ax4.plot(dv, color='red')
    ax4.set_xlim([0, 422])
    ax4.set_ylabel('DVARS')
    ax4.set_xticklabels([])
    ax4.yaxis.set_tick_params(labelsize=5)
    ax4.grid(True)

    #plot detrend
    n1  = timeseries(residuals_dt, func_gm)
    ax5 = plt.subplot2grid((8,2), (3,0),  colspan = 2, rowspan =1)
    ax5.imshow(n1, interpolation = 'none', aspect = 'auto', cmap=cm.gray, vmin =-50, vmax = 50)
    ax5.set_title('Detrend', fontsize = 8 )
    ax5.axes.get_xaxis().set_visible(False)
    ax5.axes.get_yaxis().set_visible(False)

    #plot compcor
    n2  = timeseries(residuals_cc , func_gm)
    ax6 = plt.subplot2grid((8,2), (4,0),  colspan = 2, rowspan =1)
    ax6.imshow(n2, interpolation = 'none', aspect = 'auto', cmap=cm.gray, vmin =-50, vmax = 50)
    ax6.set_title('Detrend + Compocor + Friston 24 + BP', fontsize = 8 )
    ax6.axes.get_xaxis().set_visible(False)
    ax6.axes.get_yaxis().set_visible(False)

    # #plot global
    n3  = timeseries(residuals_gl, func_gm)
    ax7 = plt.subplot2grid((8,2), (5,0),  colspan = 2, rowspan =1)
    ax7.imshow(n3, interpolation = 'none', aspect = 'auto', cmap=cm.gray, vmin =-50, vmax = 50)
    ax7.set_title('Detrend + Global + WMCSF + Friston 24 + BP', fontsize = 8)
    ax7.axes.get_xaxis().set_visible(False)
    ax7.axes.get_yaxis().set_visible(False)

    # aroma cc
    n4  = timeseries(aroma_cc, func_gm)
    ax7 = plt.subplot2grid((8,2), (6,0),  colspan = 2, rowspan =1)
    ax7.imshow(n4, interpolation = 'none', aspect = 'auto', cmap=cm.gray, vmin =-50, vmax = 50)
    ax7.set_title('AROMA + Detrend + Compocor + Friston 24 + BP', fontsize = 8)
    ax7.axes.get_xaxis().set_visible(False)
    ax7.axes.get_yaxis().set_visible(False)

    # aroma gobal
    n5  = timeseries(aroma_gl, func_gm)
    ax7 = plt.subplot2grid((8,2), (7,0),  colspan = 2, rowspan =1)
    ax7.imshow(n5, interpolation = 'none', aspect = 'auto', cmap=cm.gray, vmin =-50, vmax = 50)
    ax7.set_title('AROMA + Detrend + Global + WMCSF + Friston 24 + BP', fontsize = 8)
    ax7.axes.get_xaxis().set_visible(False)
    ax7.axes.get_yaxis().set_visible(False)

    png_name = str(os.path.join(os.getcwd(), out_name))
    plt.savefig(png_name, dpi=190, bbox_inches='tight')
    plt.close()

    return fig

def _calc_rows_columns(ratio, n_images):
    rows = 1
    for _ in range(100):
        columns = math.floor(ratio * rows)
        total = rows * columns
        if total > n_images:
            break

        columns = math.ceil(ratio * rows)
        total = rows * columns
        if total > n_images:
            break
        rows += 1
    return rows, columns

def plot_mosaic(nifti_file, output_name , title=None, overlay_mask = None, figsize=(11.7,8.3)):
    if isinstance(nifti_file,str):
        nii = nb.load(nifti_file)
        mean_data = nii.get_data()
    else:
        mean_data = nifti_file

    n_images = mean_data.shape[2]
    row, col = _calc_rows_columns(figsize[0]/figsize[1], n_images)

    if overlay_mask:
        overlay_data = nb.load(overlay_mask).get_data()

    # create figures
    fig = Figure(figsize=figsize)
    FigureCanvas(fig)

    fig.subplots_adjust(top=0.85)
    for image in (range(n_images)):
        ax = fig.add_subplot(row, col, image+1)
        data_mask = np.logical_not(np.isnan(mean_data))
        if overlay_mask:
            ax.set_rasterized(True)
        ax.imshow(np.fliplr(mean_data[:,:,image].T), vmin=np.percentile(mean_data[data_mask], 0.5),
                   vmax=np.percentile(mean_data[data_mask],99.5),
                   cmap=cm.Greys_r, interpolation='nearest', origin='lower')  # @UndefinedVariable
        if overlay_mask:
            cmap = cm.Reds  # @UndefinedVariable
            cmap._init()
            alphas = np.linspace(0, 0.75, cmap.N+3)
            cmap._lut[:,-1] = alphas
            ax.imshow(np.fliplr(overlay_data[:,:,image].T), vmin=0, vmax=1,
                   cmap=cmap, interpolation='nearest', origin='lower')  # @UndefinedVariable

        ax.axis('off')
    fig.subplots_adjust(left = 0.05, right = 0.95, bottom = 0.05, top = 0.95, wspace=0.01, hspace=0.1)

    if not title:
        _, title = os.path.split(nifti_file)
        title += " (last modified: %s)"%time.ctime(os.path.getmtime(nifti_file))
    fig.suptitle(title, fontsize='10')

    fig.savefig(output_name)
    return fig

def _get_values_inside_a_mask(main_file, mask_file):
    main_nii = nb.load(main_file)
    main_data = main_nii.get_data()
    nan_mask = np.logical_not(np.isnan(main_data))
    mask = nb.load(mask_file).get_data() > 0

    data = main_data[np.logical_and(nan_mask, mask)]
    return data

def get_median_distribution(main_files, mask_files):
    medians = []
    for main_file, mask_file in zip(main_files, mask_files):
        med = np.median(_get_values_inside_a_mask(main_file, mask_file))
        medians.append(med)
    return medians

def plot_distrbution_of_values(main_file, mask_file, xlabel, outname, distribution=None, xlabel2=None, figsize=(11.7,8.3)):

    data = _get_values_inside_a_mask(main_file, mask_file)

    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.3)
    fig.set_size_inches(12, 8)

    gs = GridSpec(2, 1)
    ax = plt.subplot(gs[0, 0])
    sns.distplot(data.astype(np.double), kde=False, bins=100, ax=ax)
    ax.set_xlabel(xlabel)

    ax = plt.subplot(gs[1, 0])
    sns.distplot(np.array(distribution).astype(np.double), ax=ax)
    cur_val = np.median(data)
    label = "%g"%cur_val
    plot_vline(cur_val, label, ax=ax)
    ax.set_xlabel(xlabel2)

    png_name = str(os.path.join(os.getcwd(), outname))
    plt.savefig(png_name, dpi=190, bbox_inches='tight')

    return fig


def find_cut_coords(img, mask=None, activation_threshold=None):
    import warnings
    import numpy as np
    from scipy import ndimage
    from nilearn._utils import as_ndarray, new_img_like
    from nilearn._utils.ndimage import largest_connected_component
    from nilearn._utils.extmath import fast_abs_percentile
    """ Find the center of the largest activation connected component.
        Parameters
        -----------
        img : 3D Nifti1Image
            The brain map.
        mask : 3D ndarray, boolean, optional
            An optional brain mask.
        activation_threshold : float, optional
            The lower threshold to the positive activation. If None, the
            activation threshold is computed using the 80% percentile of
            the absolute value of the map.
        Returns
        -------
        x : float
            the x world coordinate.
        y : float
            the y world coordinate.
        z : float
            the z world coordinate.
    """
    data = img.get_data()
    # To speed up computations, we work with partial views of the array,
    # and keep track of the offset
    offset = np.zeros(3)

    # Deal with masked arrays:
    if hasattr(data, 'mask'):
        not_mask = np.logical_not(data.mask)
        if mask is None:
            mask = not_mask
        else:
            mask *= not_mask
        data = np.asarray(data)

    # Get rid of potential memmapping
    data = as_ndarray(data)
    my_map = data.copy()
    if mask is not None:
        slice_x, slice_y, slice_z = ndimage.find_objects(mask)[0]
        my_map = my_map[slice_x, slice_y, slice_z]
        mask = mask[slice_x, slice_y, slice_z]
        my_map *= mask
        offset += [slice_x.start, slice_y.start, slice_z.start]

    # Testing min and max is faster than np.all(my_map == 0)
    if (my_map.max() == 0) and (my_map.min() == 0):
        return .5 * np.array(data.shape)
    if activation_threshold is None:
        activation_threshold = fast_abs_percentile(my_map[my_map != 0].ravel(),
                                                   80)
    mask = np.abs(my_map) > activation_threshold - 1.e-15
    # mask may be zero everywhere in rare cases
    if mask.max() == 0:
        return .5 * np.array(data.shape)
    mask = largest_connected_component(mask)
    slice_x, slice_y, slice_z = ndimage.find_objects(mask)[0]
    my_map = my_map[slice_x, slice_y, slice_z]
    mask = mask[slice_x, slice_y, slice_z]
    my_map *= mask
    offset += [slice_x.start, slice_y.start, slice_z.start]

    # For the second threshold, we use a mean, as it is much faster,
    # althought it is less robust
    second_threshold = np.abs(np.mean(my_map[mask]))
    second_mask = (np.abs(my_map) > second_threshold)
    if second_mask.sum() > 50:
        my_map *= largest_connected_component(second_mask)
    cut_coords = ndimage.center_of_mass(np.abs(my_map))
    x_map, y_map, z_map = cut_coords + offset

    coords = []
    coords.append(x_map)
    coords.append(y_map)
    coords.append(z_map)

    # Return as a list of scalars
    return coords

def plot_3d_overlay(underlay_file, overlay_file, out_filename, dpi):
    import nibabel as nb

    import matplotlib


    underlay = nb.load(underlay_file).get_data()
    overlay = nb.load(overlay_file).get_data()


    coords = find_cut_coords(nb.load(overlay_file))

    # convert zeros to nans for visualization purposes
    overlay[overlay==0]=np.nan

    # plot voxel on anat
    fig =plt.figure()
    fig.set_size_inches(6.5, 6.5)
    fig.subplots_adjust(wspace=0.005)
    import mpl_toolkits.axisartist.floating_axes as floating_axes
    from matplotlib.transforms import Affine2D

    #1
    ax1 = plt.subplot2grid((3,3), (0,0),  colspan = 1, rowspan =1)
    ax1.imshow(np.rot90(underlay[:,coords[1]-10,:], 3), matplotlib.cm.gray)
    ax1.imshow(np.rot90(overlay[:,coords[1]-10,:],3 ) , matplotlib.cm.rainbow_r, alpha = 0.5, origin='lower')
    ax1.set_xlim(23, 65)
    ax1.set_ylim(15,50)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)

    #2
    ax3 = plt.subplot2grid((3,3), (0,1),  colspan = 1, rowspan =1)
    ax3.imshow(np.rot90(underlay[:,coords[1]-7,:], 3), matplotlib.cm.gray)
    ax3.imshow(np.rot90(overlay[:,coords[1]-7,:],3 ) , matplotlib.cm.rainbow_r, alpha = 0.5, origin='lower')
    ax3.set_xlim(23, 65)
    ax3.set_ylim(15,50)
    ax3.axes.get_yaxis().set_visible(False)
    ax3.axes.get_xaxis().set_visible(False)


    #2
    ax3 = plt.subplot2grid((3,3), (0,2),  colspan = 1, rowspan =1)
    ax3.imshow(np.rot90(underlay[:,coords[1]-3,:], 3), matplotlib.cm.gray)
    ax3.imshow(np.rot90(overlay[:,coords[1]-3,:],3 ) , matplotlib.cm.rainbow_r, alpha = 0.5, origin='lower')
    ax3.set_xlim(23, 65)
    ax3.set_ylim(15,50)
    ax3.axes.get_yaxis().set_visible(False)
    ax3.axes.get_xaxis().set_visible(False)


    #3
    ax4 = plt.subplot2grid((3,3), (1,0),  colspan = 1, rowspan =1)
    ax4.imshow(np.rot90(underlay[:,coords[1],:], 3), matplotlib.cm.gray)
    ax4.imshow(np.rot90(overlay[:,coords[1],:],3 ) , matplotlib.cm.rainbow_r, alpha = 0.5, origin='lower')
    ax4.set_xlim(23, 65)
    ax4.set_ylim(15,50)
    ax4.axes.get_yaxis().set_visible(False)
    ax4.axes.get_xaxis().set_visible(False)

    #3
    ax5 = plt.subplot2grid((3,3), (1,1),  colspan = 1, rowspan =1)
    ax5.imshow(np.rot90(underlay[:,coords[1]+3,:], 3), matplotlib.cm.gray)
    ax5.imshow(np.rot90(overlay[:,coords[1]+3,:],3 ) , matplotlib.cm.rainbow_r, alpha = 0.4, origin='lower')
    ax5.set_xlim(23, 65)
    ax5.set_ylim(15,50)
    ax5.axes.get_yaxis().set_visible(False)
    ax5.axes.get_xaxis().set_visible(False)

    #3
    ax6 = plt.subplot2grid((3,3), (1,2),  colspan = 1, rowspan =1)
    ax6.imshow(np.rot90(underlay[:,coords[1]+5,:], 3), matplotlib.cm.gray)
    ax6.imshow(np.rot90(overlay[:,coords[1]+5,:],3 ) , matplotlib.cm.rainbow_r, alpha = 0.4, origin='lower')
    ax6.set_xlim(23, 65)
    ax6.set_ylim(15,50)
    ax6.axes.get_yaxis().set_visible(False)
    ax6.axes.get_xaxis().set_visible(False)

    #3
    ax7 = plt.subplot2grid((3,3), (2,0),  colspan = 1, rowspan =1)
    ax7.imshow(np.rot90(underlay[:,coords[1]+7,:], 3), matplotlib.cm.gray)
    ax7.imshow(np.rot90(overlay[:,coords[1]+7,:],3 ) , matplotlib.cm.rainbow_r, alpha = 0.4, origin='lower')
    ax7.set_xlim(23, 65)
    ax7.set_ylim(15,50)
    ax7.axes.get_yaxis().set_visible(False)
    ax7.axes.get_xaxis().set_visible(False)

    #3
    ax8 = plt.subplot2grid((3,3), (2,1),  colspan = 1, rowspan =1)
    ax8.imshow(np.rot90(underlay[:,coords[1]+10,:], 3), matplotlib.cm.gray)
    ax8.imshow(np.rot90(overlay[:,coords[1]+10,:],3 ) , matplotlib.cm.rainbow_r, alpha = 0.4, origin='lower')
    ax8.set_xlim(23, 65)
    ax8.set_ylim(15,50)
    ax8.axes.get_yaxis().set_visible(False)
    ax8.axes.get_xaxis().set_visible(False)

    #3
    ax9 = plt.subplot2grid((3,3), (2,2),  colspan = 1, rowspan =1)
    ax9.imshow(np.rot90(underlay[:,coords[1]+12,:], 3), matplotlib.cm.gray)
    ax9.imshow(np.rot90(overlay[:,coords[1]+12,:],3 ) , matplotlib.cm.rainbow_r, alpha = 0.4, origin='lower')
    ax9.set_xlim(23, 65)
    ax9.set_ylim(15,50)
    ax9.axes.get_yaxis().set_visible(False)
    ax9.axes.get_xaxis().set_visible(False)


    fig.tight_layout()


    fig.savefig(out_filename, dpi=dpi, bbox_inches='tight')
