from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import lsst.afw.math as afw_math
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from lsst.eo.pipe.plotting import plot_focal_plane


__all__ = ['apply_isr', 'compute_ctis', 'EperTask', 'EperFpPlotsTask']


def get_overscan_region(amp_info, dx):
    """
    Get the bounding box of the serial overscan region omitting dx
    leading and trailing edge pixels.
    """
    bbox = amp_info.getRawHorizontalOverscanBBox()
    minY = amp_info.getRawBBox().minY
    height = amp_info.getRawBBox().height
    return lsst.geom.Box2I(lsst.geom.Point2I(bbox.minX + dx, minY),
                           lsst.geom.Extent2I(bbox.width - 2*dx, height))


def apply_overscan_correction(raw, amp, dx=4, method='1d_poly', deg=2):
    """
    Apply serial overcan correction to the specified amp.
    This function will modify the raw pixel data in place.

    It returns the full segment, including pre- and overscan regions,
    after applying the correction to all rows.
    """
    raw_image = raw.getImage()
    det = raw.getDetector()
    amp_info = det[amp]
    full_segment = raw_image.Factory(raw_image, amp_info.getRawBBox())
    overscan = raw_image.Factory(raw_image, get_overscan_region(amp_info, dx))
    if method == 'median':
        full_segment -= np.median(overscan.array)
    elif method == 'mean':
        full_segment -= np.mean(overscan.array)
    elif method in ('median_by_row', '1d_poly'):
        row_vals = np.median(overscan.array, axis=1)
        ny = full_segment.array.shape[0]
        if method == '1d_poly':
            rows = np.arange(ny)
            func = np.poly1d(np.polyfit(rows, row_vals, deg))
            row_vals = func(rows)
        for j in range(ny):
            full_segment.array[j, :] -= row_vals[j]
    return full_segment


def subtract_calib(raw, calib, amp, scale=1):
    """
    Subtract the calibration image for the specifed amp.  The scale option
    is used to multiply the calib image, i.e., if it is a dark and scale is
    the dark time.

    This function modifies the raw pixel data in place.
    """
    raw_image = raw.getImage()
    amp_info = raw.getDetector()[amp]
    calib_image = calib.getImage()
    if scale != 1:
        calib_image = calib_image.Factory(calib_image, deep=True)
        calib_image *= scale
    calib_amp_info = calib.getDetector()[amp]
    raw_image[amp_info.getRawDataBBox()] \
        -= calib_image[calib_amp_info.getBBox()]


def apply_isr(raw, bias, dark, amp, dx=4, oscan_method='median_by_row', deg=2):
    """
    Apply a bare-bones ISR consisting of serial overscan correction,
    combined bias subtraction, combined dark correction.

    This function modifies the raw pixel data in place.

    Return the full amplifier segment.
    """
    full_segment = apply_overscan_correction(raw, amp, dx=4,
                                             method=oscan_method, deg=deg)
    subtract_calib(raw, bias, amp)
    subtract_calib(raw, dark, amp, scale=raw.getMetadata().get('DARKTIME'))
    return full_segment


def compute_ctis(processed_segment, raw_amp_info, npix=3):
    """
    For an ISR-processed amplifier segment, compute the serial and
    parallel CTI using the EPER method.
    """
    imarr = processed_segment.array
    if raw_amp_info.getRawFlipX():
        imarr[:] = imarr[:, ::-1]
    if raw_amp_info.getRawFlipY():
        imarr[:] = imarr[::-1, :]
    bbox = raw_amp_info.getRawDataBBox()
    # Serial CTI
    ncol = bbox.width
    lastcol = bbox.maxX
    signal = np.sum(imarr[:, lastcol])
    trailed = np.sum(imarr[:, lastcol + 1:lastcol + 1 + npix])
    scti = (trailed/signal)/(ncol + 1)
    # Parallel CTI
    nrow = bbox.height
    lastrow = bbox.maxY
    signal = np.sum(imarr[lastrow, :])
    trailed = np.sum(imarr[lastrow + 1:lastrow + 1 + npix, :])
    pcti = (trailed/signal)/(nrow + 1)
    return scti, pcti


class EperTaskConnections(pipeBase.PipelineTaskConnections,
                          dimensions=("instrument", "detector")):
    raws = cT.Input(
        name="raw",
        doc="Raw pixel data taken for combined flat generation",
        storageClass="Exposure",
        dimensions=("instrument", "detector", "exposure"),
        multiple=True,
        deferLoad=True)

    bias = cT.Input(
        name="bias_frame",
        doc="Combined bias frame",
        storageClass="Exposure",
        dimensions=("instrument", "detector"),
        isCalibration=True)

    dark = cT.Input(
        name="bias_frame",
        doc="Combined bias frame",
        storageClass="Exposure",
        dimensions=("instrument", "detector"),
        isCalibration=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    eper_stats = cT.Output(
        name="eper_stats",
        doc="Serial and parallel CTI values for each amp in a CCD.",
        storageClass="DataFrame",
        dimensions=("instrument", "detector"))


class EperTaskConfig(pipeBase.PipelineTaskConfig,
                     pipelineConnections=EperTaskConnections):
    nx_skip = pexConfig.Field(
        doc=("Number columns at leading and trailing edge of "
             "the serial overscan to omit when estimating the "
             "serial overscan correction."),
        default=4,
        dtype=int)
    overscan_pixels = pexConfig.Field(
        doc=("Number of overscan rows or columns to use for "
             "evaluating the trailed signal in the overscan regions."),
        default=3,
        dtype=int)
    oscan_method = pexConfig.ChoiceField(
        doc="Overscan modeling method",
        default="1d_poly",
        dtype=str,
        allowed={
            "mean": "Mean of all selected pixels in overscan region",
            "median": "Median of all selected pixels in overscan region",
            "median_by_row": "Median of each row of selected pixels.",
            "1d_poly": "1D polynomial of degree 2 fit to median_by_row data"})
    polynomial_degree = pexConfig.Field(
        doc="Degree of polynomial to fit to overscan row medians",
        default=2,
        dtype=int)
    max_raws = pexConfig.Field(
        doc="Maximum number of raw flats to included in the combined flat.",
        default=10,
        dtype=int)


class EperTask(pipeBase.PipelineTask):

    ConfigClass = EperTaskConfig
    _DefaultName = "eperTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.dx = self.config.nx_skip
        self.npix = self.config.overscan_pixels
        self.oscan_method = self.config.oscan_method
        self.deg = self.config.polynomial_degree
        self.max_raws = self.config.max_raws

    def run(self, raws, bias, dark, camera):
        det = camera[raws[0].dataId['detector']]
        det_name = det.getName()
        data = defaultdict(list)
        for amp, amp_info in enumerate(det):
            amp_name = amp_info.getName()
            images = []
            for handle in raws[:self.max_raws]:
                raw = handle.get()
                images.append(apply_isr(raw, bias, dark, amp, dx=self.dx,
                                        oscan_method=self.oscan_method,
                                        deg=self.deg))
            combined_flat = afw_math.statisticsStack(images, afw_math.MEDIAN)
            scti, pcti = compute_ctis(combined_flat, det[amp], npix=self.npix)
            data['det_name'].append(det_name)
            data['amp_name'].append(amp_name)
            data['scti'].append(scti)
            data['pcti'].append(pcti)
        return pipeBase.Struct(eper_stats=pd.DataFrame(data))


class EperFpPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                 dimensions=("instrument",)):
    eper_stats = cT.Input(
        name="eper_stats",
        doc="Serial and parallel CTI values obtained using EPER method",
        storageClass="DataFrame",
        dimensions=("instrument", "detector"),
        multiple=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    scti_eper_plot = cT.Output(
        name="scti_eper_plot",
        doc="Focal plane mosaic of serial CTI for each amp",
        storageClass="Plot",
        dimensions=("instrument",))

    pcti_eper_plot = cT.Output(
        name="[cti_eper_plot",
        doc="Focal plane mosaic of parallel CTI for each amp",
        storageClass="Plot",
        dimensions=("instrument",))


class EperFpPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                            pipelineConnections=EperFpPlotsTaskConnections):
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    zmin = pexConfig.Field(doc="Minimum of color bar range.",
                           dtype=float, default=0)
    zmax = pexConfig.Field(doc="Maximum of color bar range.",
                           dtype=float, default=1e-5)
    zscale_factor = pexConfig.Field(doc=("Scale factor to apply to z-values."
                                         "This should be a float-convertable "
                                         "string so that formatting is taken "
                                         "care of automatically when "
                                         "rendering the label in matplotlib"),
                                    dtype=str, default="1e-6")


class EperFpPlotsTask(pipeBase.PipelineTask):
    """
    Create a focal plane mosaic of the serial and parallel CTI results.
    """
    ConfigClass = EperFpPlotsTaskConfig
    _defaultName = 'eperFpPlotsTask'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.yfigsize, self.config.xfigsize
        self.z_range = self.config.zmin, self.config.zmax
        self.zscale_factor = self.config.zscale_factor

    def run(self, eper_stats, camera):
        # Unpack the divisadero data for plotting.
        scti_data = defaultdict(dict)
        pcti_data = defaultdict(dict)
        for handle in eper_stats:
            det = camera[handle.dataId['detector']]
            det_name = det.getName()
            df = handle.get()
            for _, row in df.iterrows():
                amp_name = row['amp_name']
                scti_data[det_name][amp_name] = row['scti']
                pcti_data[det_name][amp_name] = row['pcti']

        scti_eper_plot = plt.figure(figsize=self.figsize)
        ax = scti_eper_plot.add_subplot(111)
        plot_focal_plane(ax, scti_data, camera=camera, z_range=self.z_range,
                         scale_factor=self.zscale_factor,
                         title='Serial CTI', )

        pcti_eper_plot = plt.figure(figsize=self.figsize)
        ax = pcti_eper_plot.add_subplot(111)
        plot_focal_plane(ax, pcti_data, camera=camera, z_range=self.z_range,
                         scale_factor=self.zscale_factor,
                         title='Parallel CTI', )

        return pipeBase.Struct(scti_eper_plot=scti_eper_plot,
                               pcti_eper_plot=pcti_eper_plot)
