from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from lsst.eo.pipe.plotting import plot_focal_plane


__all__ = ['DarkCurrentTask']


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

    dark_current_plot = cT.Output(
        name="dark_current_plot",
        doc="Focal plane mosaic of per-amp dark current values",
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
        amp_data = defaultdict(dict)
        data = defaultdict(list)
        for handle in darks:
            dark = handle.get()
            det = dark.getDetector()
            det_name = det.getName()
            for amp, amp_info in enumerate(det):
                amp_name = amp_info.getName()
                amp_image = dark.getImage()[amp_info.getBBox()]
                dark_current = np.percentile(amp_image.array, self.percentile)
                data['det_name'].append(det_name)
                data['amp_name'].append(amp_name)
                data['dark_current'].append(dark_current)
                amp_data[det_name][amp_name] = dark_current
        dark_current_stats = pd.DataFrame(data)
        dark_current_plot = plt.figure(figsize=self.figsize)
        ax = dark_current_plot.add_subplot(111)
        title = f"Dark current, {self.percentile} percentile [e-/s]"
        plot_focal_plane(ax, amp_data, camera=camera, z_range=self.z_range,
                         scale_factor=self.zscale_factor, title=title)
        return pipeBase.Struct(dark_current_stats=dark_current_stats,
                               dark_current_plot=dark_current_plot)