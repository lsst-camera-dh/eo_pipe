from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.daf.butler as daf_butler
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from .plotting import plot_focal_plane, hist_amp_data, append_acq_run
from .dsref_utils import get_plot_locations_by_dstype


__all__ = ['DarkCurrentTask']


def get_amp_data(repo, collections):
    """Return amp-level dark current data."""
    butler = daf_butler.Butler(repo, collections=collections)
    dsrefs = list(set(butler.registry.queryDatasets('dark_current_stats',
                                                    findFirst=True)))
    df = butler.getDirect(dsrefs[0])
    amp_data = {column: defaultdict(dict) for column in df.columns
                if column.startswith('dark_current')}
    for _, row in df.iterrows():
        for column in amp_data:
            amp_data[column][row.det_name][row.amp_name] = row[column]
    return {column: dict(data) for column, data in amp_data.items()}


def get_plot_locations(repo, collections):
    dstypes = ('dark_current_percentile_plot', 'dark_current_median_plot',
               'dark_current_percentile_hist', 'dark_current_median_hist')
    return get_plot_locations_by_dstype(repo, collections, dstypes)


class DarkCurrentTaskConnections(pipeBase.PipelineTaskConnections,
                                 dimensions=("instrument",)):
    darks = cT.Input(
        name="dark",
        doc="Combined dark frames",
        storageClass="Exposure",
        dimensions=("instrument", "detector"),
        isCalibration=True,
        multiple=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    dark_current_stats = cT.Output(
        name="dark_current_stats",
        doc="Dark current values for each amp in the focal plane",
        storageClass="DataFrame",
        dimensions=("instrument",))

    dark_current_percentile_plot = cT.Output(
        name="dark_current_percentile_plot",
        doc="Focal plane mosaic of per-amp percentile dark current values",
        storageClass="Plot",
        dimensions=("instrument",))

    dark_current_percentile_hist = cT.Output(
        name="dark_current_percentile_hist",
        doc="Histogram of per-amp percentile dark current values",
        storageClass="Plot",
        dimensions=("instrument",))

    dark_current_median_plot = cT.Output(
        name="dark_current_median_plot",
        doc="Focal plane mosaic of per-amp median dark current values",
        storageClass="Plot",
        dimensions=("instrument",))

    dark_current_median_hist = cT.Output(
        name="dark_current_median_hist",
        doc="Histogram of per-amp median dark current values",
        storageClass="Plot",
        dimensions=("instrument",))


class DarkCurrentTaskConfig(pipeBase.PipelineTaskConfig,
                            pipelineConnections=DarkCurrentTaskConnections):
    percentile = pexConfig.Field(
        doc="Percentile of dark current pixel values to compute",
        dtype=float,
        default=95.)
    xfigsize = pexConfig.Field(
        doc="Figure size x-direction in inches.",
        dtype=float,
        default=12)
    yfigsize = pexConfig.Field(
        doc="Figure size y-direction in inches.",
        dtype=float,
        default=12)
    zmin = pexConfig.Field(
        doc="Minimum of color bar range.",
        dtype=float,
        default=0)
    zmax = pexConfig.Field(
        doc="Maximum of color bar range.",
        dtype=float,
        default=0.3)
    zscale_factor = pexConfig.Field(
        doc=("Scale factor to apply to z-values.  This should be a "
             "float-convertable string so that formatting is taken "
             "care of automatically when rendering the label in matplotlib."),
        dtype=str,
        default="1")
    acq_run = pexConfig.Field(doc="Acquisition run number.",
                              dtype=str, default="")


class DarkCurrentTask(pipeBase.PipelineTask):
    """Task to measure serial and parallel CTI using the EPER method."""
    ConfigClass = DarkCurrentTaskConfig
    _DefaultName = "darkCurrentTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.percentile = self.config.percentile
        self.figsize = self.config.yfigsize, self.config.xfigsize
        self.z_range = self.config.zmin, self.config.zmax
        self.zscale_factor = self.config.zscale_factor

    def run(self, darks, camera):
        amp_data = {'percentile': defaultdict(dict),
                    'median': defaultdict(dict)}
        data = defaultdict(list)
        for handle in darks:
            dark = handle.get()
            det = dark.getDetector()
            det_name = det.getName()
            for amp, amp_info in enumerate(det):
                amp_name = amp_info.getName()
                amp_image = dark.getImage()[amp_info.getBBox()]
                dark_current_percentile = np.percentile(amp_image.array,
                                                        self.percentile)
                dark_current_median = np.median(amp_image.array)
                data['det_name'].append(det_name)
                data['amp_name'].append(amp_name)
                dc_name = f"dark_current_{int(self.percentile):02d}"
                data[dc_name].append(dark_current_percentile)
                data['dark_current_median'].append(dark_current_median)
                amp_data['percentile'][det_name][amp_name] = dark_current_percentile
                amp_data['median'][det_name][amp_name] = dark_current_median
        dark_current_stats = pd.DataFrame(data)
        dark_current_percentile_plot = plt.figure(figsize=self.figsize)
        ax = dark_current_percentile_plot.add_subplot(111)
        title = append_acq_run(self, f"Dark current, {self.percentile} "
                               "percentile [e-/s]")
        plot_focal_plane(ax, amp_data['percentile'], camera=camera,
                         z_range=self.z_range,
                         scale_factor=self.zscale_factor, title=title)
        dark_current_percentile_hist = plt.figure()
        hist_amp_data(amp_data['percentile'], "Dark current [e-/s]",
                      hist_range=self.z_range, scale_factor=self.zscale_factor,
                      title=title)

        dark_current_median_plot = plt.figure(figsize=self.figsize)
        ax = dark_current_median_plot.add_subplot(111)
        title = append_acq_run(self, "Dark current, median [e-/s]")
        plot_focal_plane(ax, amp_data['median'], camera=camera,
                         z_range=self.z_range,
                         scale_factor=self.zscale_factor, title=title)
        dark_current_median_hist = plt.figure()
        hist_amp_data(amp_data['median'], "Dark current [e-/s]",
                      hist_range=self.z_range, scale_factor=self.zscale_factor,
                      title=title)
        return pipeBase.Struct(dark_current_stats=dark_current_stats,
                               dark_current_percentile_plot=dark_current_percentile_plot,
                               dark_current_percentile_hist=dark_current_percentile_hist,
                               dark_current_median_plot=dark_current_median_plot,
                               dark_current_median_hist=dark_current_median_hist)
