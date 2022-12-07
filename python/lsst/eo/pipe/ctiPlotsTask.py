from collections import defaultdict
import matplotlib.pyplot as plt
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from lsst.eo.pipe.plotting import plot_focal_plane


__all__ = ["CtiFpPlotsTask"]


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
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    serial_cti = cT.Output(
        name="serial_cti",
        doc="Serial CTI per amp",
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

        figure = plt.figure(figsize=self.figsize)
        ax = plt.subplot(111)
        plot_focal_plane(ax, scti, camera=camera, z_range=self.zrange,
                         scale_factor=self.zscale_factor)

        return pipeBase.Struct(serial_cti=figure)
