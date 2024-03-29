from collections import defaultdict
import matplotlib.pyplot as plt
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from .plotting import plot_focal_plane, hist_amp_data, append_acq_run
from .dsref_utils import get_plot_locations_by_dstype


__all__ = ["CtiFpPlotsTask"]


def get_plot_locations(repo, collections):
    dstypes = ('serial_cti', 'serial_cti_hist')
    return get_plot_locations_by_dstype(repo, collections, dstypes)


class CtiFpPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                dimensions=("instrument",)):
    cti_calibs = cT.Input(
        name="cpCtiCalib",
        doc="CTI calibration outputs",
        storageClass="IsrCalib",
        dimensions=("instrument", "detector"),
        isCalibration=True,
        multiple=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",))

    serial_cti = cT.Output(
        name="serial_cti",
        doc="Serial CTI per amp",
        storageClass="Plot",
        dimensions=("instrument",))

    serial_cti_hist = cT.Output(
        name="serial_cti_hist",
        doc="Histogram of serial CTI values per amp",
        storageClass="Plot",
        dimensions=("instrument",))


class CtiFpPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                           pipelineConnections=CtiFpPlotsTaskConnections):
    """Configuration for CtiFpPlotsTask."""
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
                                         "care of automatically when rendering "
                                         "the label in matplotlib"),
                                    dtype=str, default="1e-6")


class CtiFpPlotsTask(pipeBase.PipelineTask):
    """Create focal plane mosaic of serial CTI per amp."""
    ConfigClass = CtiFpPlotsTaskConfig
    _DefaultName = "ctiFpPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.yfigsize, self.config.xfigsize
        self.zrange = self.config.zmin, self.config.zmax
        self.zscale_factor = self.config.zscale_factor

    def run(self, cti_calibs, camera):
        scti = defaultdict(dict)
        for handle in cti_calibs:
            det = camera[handle.dataId['detector']]
            det_name = det.getName()
            for amp in det:
                channel = amp.getName()
                calib = handle.get()
                try:
                    cti = calib.globalCti[channel]
                    scti[det_name][channel] = cti
                except KeyError:
                    pass

        serial_cti = plt.figure(figsize=self.figsize)
        ax = plt.subplot(111)
        title = append_acq_run(self, 'SCTI from deferredCharge task')
        plot_focal_plane(ax, scti, camera=camera, z_range=self.zrange,
                         scale_factor=self.zscale_factor, title=title)

        serial_cti_hist = plt.figure()
        hist_amp_data(scti, 'serial CTI', hist_range=self.zrange,
                      scale_factor=self.zscale_factor, title=title)

        return pipeBase.Struct(serial_cti=serial_cti,
                               serial_cti_hist=serial_cti_hist)
