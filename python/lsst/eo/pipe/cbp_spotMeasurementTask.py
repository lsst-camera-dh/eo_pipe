from lsst.afw import cameraGeom
import sys
import numpy as np
import lsst.afw.table as afw_table
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
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
        dimensions=("instrument",))

    spot_catalog = cT.Output(
        name="spot_catalog",
        doc="Catalog of spot measurements.",
        storageClass="SourceCatalog",
        dimensions=("instrument", "detector"))


class SpotMeasurementTaskConfig(pipeBase.PipelineTaskConfig,
                                pipelineConnections=SpotMeasurementTaskConnections):
    """Configuration for SpotMeasurementTask."""
    threshold_adu = pexConfig.Field(doc="Spot threshold in ADU",
                                    dtype=float, default=10.0) #10.0 default
    minarea = pexConfig.Field(doc="Minimum pixel area for spots.",
                              dtype=float, default=20000.0) # 30000 default
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
        sys.path.append('/sdf/data/rubin/user/amouroux/comissioning/cbp_codes/scripts/spot_measurement/')
        from photometry import AperturePhotometry
        #from lsst.eo.throughput import throughput

        schema = afw_table.SourceTable.makeMinimalSchema()
        schema.addField("det_name", str, doc="Detector name", size=10)
        schema.addField("exposure", "L", doc="Exposure ID")
        schema.addField("x", "F", doc="x-coordinate of spot centroid on CCD")
        schema.addField("y", "F", doc="y-coordinate of spot centroid on CCD")
        schema.addField("x_fp", "F", doc="x-coordinate of spot centroid in "
                        "focal plane coordinates (mm)")
        schema.addField("y_fp", "F", doc="y-coordinate of spot centroid in "
                        "focal plane coordinates (mm)")
        schema.addField("radius", "F", doc="width of spot in pixels")
        schema.addField("signal", "F", doc="Spot signal")
        table = afw_table.SourceTable.make(schema)
        catalog = afw_table.SourceCatalog(table)
        for handle in exposure_handles:
            det = camera[handle.dataId['detector']]
            pix_to_fp = det.getTransform(cameraGeom.PIXELS,
                                         cameraGeom.FOCAL_PLANE)
            det_name = det.getName()
            exposure = handle.get()
            exp_id = handle.dataId['exposure']
            
            ap = AperturePhotometry()
            signal = ap.do_aperture_photometry()
            x, y = ap.x, ap.y
            radius = ap.radiuses
            fpx, fpy = pix_to_fp.getMapping().applyForward(np.vstack((x, y)))
            for i in range(len(ap.x)):
                record = catalog.addNew()
                record.set('det_name', det_name)
                record.set('exposure', exp_id)
                record.set('x', x[i])
                record.set('y', y[i])
                record.set('x_fp', fpx[i])
                record.set('y_fp', fpy[i])
                record.set('radius', radius[i])
                record.set('signal', signal[i])
            #spot_info = throughput.get_spots_counts(
            #    exposure,
            #    threshold_adu=self.config.threshold_adu,
            #    minarea=self.config.minarea,
            #    maxarea=self.config.maxarea,
            #    force_circle=self.config.force_circle)

            #for signal, (x, y), fp in spot_info:
            #    fp_pos = pix_to_fp.applyForward(lsst.geom.Point2D(x, y))
            #    record = catalog.addNew()
            #    record.set('det_name', det_name)
            #    record.set('exposure', exp_id)
            #    record.set('x', x)
            #    record.set('y', y)
            #    record.set('x_fp', fp_pos.x)
            #    record.set('y_fp', fp_pos.y)
            #    record.set('signal', signal)
            #    record.setFootprint(fp)
        return pipeBase.Struct(spot_catalog=catalog)
