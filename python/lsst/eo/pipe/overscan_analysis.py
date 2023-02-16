import itertools
import numpy as np
import astropy.visualization as viz
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable


__all__ = ('raft_oscan_correlations',)


def raft_oscan_correlations(bias_files, buffer=10, title='',
                            vrange=None, stretch=viz.LinearStretch,
                            figsize=(8, 8)):
    """
    Compute the correlation coefficients between the overscan pixels
    of the 144 amplifiers in raft.

    Parameters
    ----------
    bias_files: dict
        Dictionary of bias image files, indexed by sensor slot id.
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
        raw = ref.get()
        det = raw.getDetector()
        for amp in det:
            bbox = amp.getRawSerialOverscanBBox()
            bbox.grow(-buffer)
            overscans.append(raw.getImage()[bbox].array)
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
    image = ax.imshow(data, interpolation='none', norm=norm)
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
