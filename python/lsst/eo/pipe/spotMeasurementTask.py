import lsst.afw.table as afw_table
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from lsst.eo.throughput import throughput


__all__ = ["SpotMeasurementTask"]


class SpotMeasurementTaskConnections(pipeBase.PipelineTaskConnections,
                                     dimensions=("instrument", "detector")):
    exposure_handles = cT.Input(
        name="postISRCCD",
        doc="ISR'd exposures to analyze",
        storageClass="Exposure",
        dimensions=("instrument", "exposure", "detector"),
        multiple=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    spot_catalog = cT.Output(
        name="spot_catalog",
        doc="Catalog of spot measurements.",
        storageClass="SourceCatalog",
        dimensions=("instrument", "detector"))


class SpotMeasurementTaskConfig(pipeBase.PipelineTaskConfig,
                                pipelineConnections=SpotMeasurementTaskConnections):
    """Configuration for SpotMeasurementTask."""
    threshold_adu = pexConfig.Field(doc="Spot threshold in ADU",
                                    dtype=float, default=10.0)
    minarea = pexConfig.Field(doc="Minimum pixel area for spots.",
                              dtype=float, default=30000.0)
    maxarea = pexConfig.Field(doc="Maximum pixel area for spots.",
                              dtype=float, default=4000000.0)
    force_circle = pexConfig.Field(doc="Force spot shape to be a circle",
                                   dtype=bool, default=True)


class SpotMeasurementTask(pipeBase.PipelineTask):
    """Spot measurement task for CCOB narrow bean data."""
    ConfigClass = SpotMeasurementTaskConfig
    _DefaultName = "spotMeasurementTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, exposure_handles, camera):
        schema = afw_table.SourceTable.makeMinimalSchema()
        schema.addField("det_name", str, doc="Detector name", size=10)
        schema.addField("exposure", "I", doc="Exposure ID")
        schema.addField("x", "F", doc="x-coordinate of spot")
        schema.addField("y", "F", doc="y-coordinate of spot")
        schema.addField("signal", "F", doc="Spot signal")
        table = afw_table.SourceTable.make(schema)
        catalog = afw_table.SourceCatalog(table)
        for handle in exposure_handles:
            det_name = camera[handle.dataId['detector']].getName()
            exposure = handle.get()
            exp_id = handle.dataId['exposure']
            spot_info = throughput.get_spots_counts(
                exposure,
                threshold_adu=self.config.threshold_adu,
                minarea=self.config.minarea,
                maxarea=self.config.maxarea,
                force_circle=self.config.force_circle)
            for signal, (x, y), fp in spot_info:
                record = catalog.addNew()
                record.set('det_name', det_name)
                record.set('exposure', exp_id)
                record.set('x', x)
                record.set('y', y)
                record.set('signal', signal)
                record.setFootprint(fp)
        return pipeBase.Struct(spot_catalog=catalog)
