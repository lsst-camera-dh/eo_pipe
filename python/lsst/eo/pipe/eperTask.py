from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import lsst.afw.math as afw_math
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.daf.butler as daf_butler
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from .plotting import plot_focal_plane, hist_amp_data, append_acq_run
from .isr_utils import apply_minimal_isr
from .dsref_utils import get_plot_locations_by_dstype


__all__ = ['compute_ctis', 'EperTask', 'EperFpPlotsTask']


def get_amp_data(repo, collections):
    """Get EPER results for each amp in the focal plane."""
    butler = daf_butler.Butler(repo, collections=collections)
    dsrefs = list(set(butler.registry.queryDatasets('eper_stats',
                                                    findFirst=True)))
    amp_data = defaultdict(lambda: defaultdict(dict))
    for dsref in dsrefs:
        df = butler.get(dsref)
        for _, row in df.iterrows():
            for field in ('scti', 'pcti'):
                amp_data[field][row.det_name][row.amp_name] = row[field]
    return {field: dict(data) for field, data in amp_data.items()}


def get_plot_locations(repo, collections):
    dstypes = ('scti_eper_plot', 'pcti_eper_plot',
               'scti_eper_hist', 'pcti_eper_hist')
    return get_plot_locations_by_dstype(repo, collections, dstypes)


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
    firstcol = bbox.minX
    lastcol = bbox.maxX
    firstrow = bbox.minY
    lastrow = bbox.maxY
    # Serial CTI
    ncol = bbox.width
    signal = np.sum(imarr[firstrow:lastrow + 1, lastcol])
    trailed = np.sum(imarr[firstrow:lastrow + 1,
                           lastcol + 1:lastcol + 1 + npix])
    scti = (trailed/signal)/(ncol + 1)
    # Parallel CTI
    nrow = bbox.height
    signal = np.sum(imarr[lastrow, firstcol:lastcol+1])
    trailed = np.sum(imarr[lastrow + 1:lastrow + 1 + npix, firstcol:lastcol+1])
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
        name="dark_frame",
        doc="Combined dark frame",
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
        doc=("Number columns at the leading and trailing edges of "
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
            "median_per_row": "Median of each row of selected pixels",
            "1d_poly": "1D polynomial of degree 2 fit to median_per_row data"})
    polynomial_degree = pexConfig.Field(
        doc="Degree of polynomial to fit to overscan row medians",
        default=2,
        dtype=int)
    max_raws = pexConfig.Field(
        doc="Maximum number of raw flats to included in the combined flat.",
        default=10,
        dtype=int)


class EperTask(pipeBase.PipelineTask):
    """Task to measure serial and parallel CTI using the EPER method."""
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
                images.append(apply_minimal_isr(raw, bias, dark, amp,
                                                dx=self.dx,
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

    scti_eper_hist = cT.Output(
        name="scti_eper_hist",
        doc="Histogram of serial CTI for each amp",
        storageClass="Plot",
        dimensions=("instrument",))

    pcti_eper_plot = cT.Output(
        name="pcti_eper_plot",
        doc="Focal plane mosaic of parallel CTI for each amp",
        storageClass="Plot",
        dimensions=("instrument",))

    pcti_eper_hist = cT.Output(
        name="pcti_eper_hist",
        doc="Histogram of parallel CTI for each amp",
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
    acq_run = pexConfig.Field(doc="Acquisition run number.",
                              dtype=str, default="")


class EperFpPlotsTask(pipeBase.PipelineTask):
    """
    Create focal plane mosaics of the serial and parallel CTI results.
    """
    ConfigClass = EperFpPlotsTaskConfig
    _DefaultName = 'eperFpPlotsTask'

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
        xlabel = "Serial CTI"
        title = append_acq_run(self, xlabel)
        plot_focal_plane(ax, scti_data, camera=camera, z_range=self.z_range,
                         scale_factor=self.zscale_factor, title=title)
        scti_eper_hist = plt.figure()
        hist_amp_data(scti_data, xlabel, hist_range=self.z_range,
                      scale_factor=self.zscale_factor, title=title)

        pcti_eper_plot = plt.figure(figsize=self.figsize)
        ax = pcti_eper_plot.add_subplot(111)
        xlabel = "Parallel CTI"
        title = append_acq_run(self, xlabel)
        plot_focal_plane(ax, pcti_data, camera=camera, z_range=self.z_range,
                         scale_factor=self.zscale_factor, title=title)
        pcti_eper_hist = plt.figure()
        hist_amp_data(pcti_data, xlabel, hist_range=self.z_range,
                      scale_factor=self.zscale_factor, title=title)

        return pipeBase.Struct(scti_eper_plot=scti_eper_plot,
                               scti_eper_hist=scti_eper_hist,
                               pcti_eper_plot=pcti_eper_plot,
                               pcti_eper_hist=pcti_eper_hist)
