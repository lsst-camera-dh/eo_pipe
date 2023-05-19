import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import lsst.afw.math as afwMath
from lsst.afw.cameraGeom import utils as cgu


__all__ = ['cmap_range', 'ImageSource', 'make_mosaic', 'append_acq_run',
           'nsigma_range']


def nsigma_range(data, nsigma=3):
    stats = afwMath.makeStatistics(data, afwMath.MEDIAN | afwMath.STDEVCLIP
                                   | afwMath.STDEV)
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    if not np.isfinite(stdev):
        return None
    return (median - nsigma*stdev, median + nsigma*stdev)


def append_acq_run(cls_instance, title, suffix=None):
    if cls_instance.config.acq_run != -1:
        title = f"{title}, acq. run {cls_instance.config.acq_run.strip()}"
    if suffix is not None:
        title = f"{title}, {suffix}"
    return title


def cmap_range(image_array, nsig=5):
    pixel_data = np.array(image_array, dtype=float).flatten()
    stats = afwMath.makeStatistics(pixel_data,
                                   afwMath.STDEVCLIP | afwMath.MEDIAN)
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    vmin = max(min(pixel_data), median - nsig*stdev)
    vmax = min(max(pixel_data), median + nsig*stdev)
    return vmin, vmax


class ImageSource:
    isTrimmed = True
    background = 0.0

    def __init__(self, exposure_handles):
        self.exposure_handles = exposure_handles

    def getCcdImage(self, det, imageFactory, binSize=1, *args, **kwargs):
        ccdImage = self.exposure_handles[det.getId()].get().getImage()
        ccdImage = afwMath.binImage(ccdImage, binSize)
        return afwMath.rotateImageBy90(ccdImage,
                                       det.getOrientation().getNQuarter()), det


def make_mosaic(exposure_refs, camera, binSize, figsize, cmap, nsig,
                title=None):
    """
    Make a a mosaic of exposures using the
    lsst.afw.cameraGeom.utils.showCamera function.

    Parameters
    ----------
    exposure_refs : dict
        Dictionary of dataset references to the exposures, keyed by
        detector number.
    camera : lsst.afw.cameraGeom.Camera
        The LSST Camera object. This will typically be LSSTCam.
    binSize : int
        Rebinning size.
    figsize : (float, float)
        Figure size in inches.
    cmap : matplotlib.colors.Colormap
        Color map.
    nsig : float
        Number of clipped stdevs to use around the median for plotting range.
    title : str [None]
        Plot title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    detectorNameList = [camera[detector].getName() for
                        detector in exposure_refs]
    image_source = ImageSource(exposure_refs)
    mosaic = cgu.showCamera(camera, imageSource=image_source,
                            detectorNameList=detectorNameList,
                            binSize=binSize)
    mosaic = afwMath.flipImage(mosaic, flipLR=False, flipTB=True)
    imarr = mosaic.array
    my_plot = plt.figure(figsize=figsize)
    ax = my_plot.add_subplot(111)
    image = plt.imshow(imarr, interpolation='nearest', cmap=cmap)
    vmin, vmax = cmap_range(imarr, nsig=nsig)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    image.set_norm(norm)
    plt.tick_params(axis='both', which='both', top=False,
                    bottom=False, left=False, right=False,
                    labelbottom=False, labelleft=False)
    if title is not None:
        plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(image, cax=cax)
    return my_plot
