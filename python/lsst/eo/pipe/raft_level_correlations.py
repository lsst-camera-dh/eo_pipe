import itertools
import numpy as np
import astropy.visualization as viz
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from .isr_utils import apply_minimal_isr


__all__ = ('raft_oscan_correlations', 'raft_imaging_correlations')


def pearson_r(x, y):
    """
    Compute Pearson correlation coefficient between two arrays.
    """
    return (np.mean(x*y) - np.mean(x)*np.mean(y))/np.std(x)/np.std(y)


def raft_oscan_correlations(bias_refs, camera, buffer=10, title='',
                            vrange=None, stretch=viz.LinearStretch,
                            cmap='jet', figsize=(8, 8)):
    """
    Compute the correlation coefficients between the overscan pixels
    of the 144 amplifiers in raft.

    Parameters
    ----------
    bias_refs: dict
        Dictionary of dataset references to bias image files, indexed
        by sensor slot id.
    camera: lsst.afw.cameraGeom.Camera
        Camera to use to find amplifier flips to put the overscan
        arrays in readout order.
    buffer: int [10]
        Buffer region around perimeter of serial overscan region to
        avoid when computing the correlation coefficients.
    title: str ['']
        Plot title.
    vrange: (float, float) [None]
        Minimum and maximum values for color scale range. If None, then
        the range of the central 98th percentile of the absolute value
        of the data is used.
    stretch: astropy.visualization.BaseStretch [LinearStretch]
        Stretch to use for the color scale.
    cmap: str ['jet']
        Matplotlib color map to use.

    Returns
    -------
    (matplotlib.figure.Figure, np.array): The figure containing the plot and
        the numpy array containing the correlation coefficients.
    """
    # Create a list of overscan np.arrays, ordered by slot and amp name
    # for all of the CCDs in a raft.
    #
    # List of slots to provide specific ordering of overscan data.
    slots = 'S00 S01 S02 S10 S11 S12 S20 S21 S22'.split()
    overscans = []

    for slot in slots:
        ref = bias_refs[slot]
        exp = ref.get()
        det = exp.getDetector()
        raw_det = camera[det.getName()]
        for amp, raw_amp in zip(det, raw_det):
            bbox = amp.getRawSerialOverscanBBox()
            bbox.grow(-buffer)
            data = exp.getImage()[bbox].array.copy()
            # Put the overscan data in readout order.
            if raw_amp.getRawFlipX():
                data = data[:, ::-1]
            if raw_amp.getRawFlipY():
                data = data[::-1, :]
            overscans.append(data)
    namps = len(overscans)
    data = np.array([np.corrcoef(overscans[i[0]].ravel(),
                                 overscans[i[1]].ravel())[0, 1]
                     for i in itertools.product(range(namps), range(namps))])
    data = data.reshape((namps, namps))
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize='medium')

    interval = viz.PercentileInterval(98.)
    if vrange is None:
        vrange = interval.get_limits(np.abs(data.ravel()))
    norm = ImageNormalize(vmin=vrange[0], vmax=vrange[1], stretch=stretch())
    image = ax.imshow(data, interpolation='none', norm=norm, cmap=cmap)
    set_ticks(ax, slots, amps=16)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(image, cax=cax)

    return fig, data


def set_ticks(ax, slots, amps=16):
    """Set the tick labels, centering the slot names between amps 1 and 16."""
    major_locs = [i*amps - 0.5 for i in range(len(slots) + 1)]
    minor_locs = [amps//2 + i*amps for i in range(len(slots))]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_tick_params(which='minor', length=0)
        axis.set_major_locator(ticker.FixedLocator(major_locs))
        axis.set_major_formatter(ticker.FixedFormatter(['']*len(major_locs)))
        axis.set_minor_locator(ticker.FixedLocator(minor_locs))
        axis.set_minor_formatter(ticker.FixedFormatter(slots))


def raft_imaging_correlations(flat1_refs, flat2_refs, camera,
                              buffer=10, title='', vrange=None,
                              stretch=viz.LinearStretch, cmap='jet',
                              figsize=(8, 8)):

    """
    Compute the correlation coefficients between the imaging section
    pixels for the difference images from a flat pair for the 144
    amplifiers in raft.

    Parameters
    ----------
    flat1_refs: dict
        Dictionary of flat1 exposure references, indexed by sensor slot id.
        These should be from the same flat pair frame as the flat2_refs.
    flat2_refs: dict
        Dictionary of flat2 exposure references, indexed by sensor slot id.
    camera: lsst.afw.cameraGeom.Camera
        Camera to use to find amplifier flips to put the overscan
        arrays in readout order.
    buffer: int [10]
        Buffer region around perimeter of imaging region to
        avoid when computing the correlation coefficients.
    title: str ['']
        Plot title.
    vrange: (float, float) [None]
        Minimum and maximum values for color scale range. If None, then
        the range of the central 98th percentile of the absolute value
        of the data is used.
    stretch: astropy.visualization.BaseStretch [LinearStretch]
        Stretch to use for the color scale.
    cmap: str ['jet']
        Matplotlib color map to use.

    Returns
    -------
    (matplotlib.figure.Figure, np.array): The figure containing the plot and
        the numpy array containing the correlation coefficients.
    """
    slots = 'S00 S01 S02 S10 S11 S12 S20 S21 S22'.split()
    segments = []

    for slot in slots:
        exp1 = flat1_refs[slot].get()
        exp2 = flat2_refs[slot].get()
        det = exp1.getDetector()
        raw_det = camera[det.getName()]
        imarrs = diff_image_arrays(exp1, exp2, buffer=buffer)
        for amp_name, data in imarrs.items():
            raw_amp = raw_det[amp_name]
            if raw_amp.getRawFlipX():
                data = data[:, ::-1]
            if raw_amp.getRawFlipY():
                data = data[::-1, :]
            segments.append(data)
    namps = len(segments)
    data = np.array([pearson_r(segments[i[0]].ravel(), segments[i[1]].ravel())
                     for i in itertools.product(range(namps), range(namps))])
    data = data.reshape((namps, namps))
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize='medium')

    interval = viz.PercentileInterval(98.)
    if vrange is None:
        vrange = interval.get_limits(np.abs(data.ravel()))
    norm = ImageNormalize(vmin=vrange[0], vmax=vrange[1], stretch=stretch())
    image = ax.imshow(data, interpolation='none', norm=norm, cmap=cmap)
    set_ticks(ax, slots, amps=16)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(image, cax=cax)

    return fig, data


def diff_image_arrays(exp1, exp2, bias=None, dark=None, buffer=10):
    """
    Compute the difference images for each amp, using the Astier
    weighting scheme to account for somewhat different exposure times,
    and return a dict of the image arrays, keyed by amp.
    """
    image_arrays = {}
    det = exp1.getDetector()
    for amp in det:
        amp_name = amp.getName()
        bbox = amp.getRawDataBBox()  # bounding box of imaging section
        bbox.grow(-buffer)
        apply_minimal_isr(exp1, bias, dark, amp_name, do_parallel=True)
        apply_minimal_isr(exp2, bias, dark, amp_name, do_parallel=True)
        imarr1 = exp1.getImage()[bbox].array.copy()
        imarr2 = exp2.getImage()[bbox].array.copy()
        median1 = np.median(imarr1)
        median2 = np.median(imarr2)
        fmedian = (median1 + median2)/2.
        imarr1 *= median2/fmedian
        imarr2 *= median1/fmedian
        image_arrays[amp_name] = imarr1 - imarr2
    return image_arrays
