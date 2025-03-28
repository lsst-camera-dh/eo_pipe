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
        name="postISRCCD",#postISRCCD
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

    def run(self, exposure_handles, camera, aperture, forced = False, ref_spots=None, repo="embargo_new"):
        sys.path.append("/sdf/home/a/amouroux/DATA/eo_pipe/python/lsst/eo/pipe")
        from photometry import AperturePhotometry, ImageData, Spot
        #from lsst.eo.throughput import throughput

        schema = afw_table.SourceTable.makeMinimalSchema()
        schema.addField("detector", str, doc="Detector name", size=10)
        schema.addField("exposure", "L", doc="Exposure ID")
        schema.addField("x", "F", doc="x-coordinate of spot centroid on CCD")
        schema.addField("y", "F", doc="y-coordinate of spot centroid on CCD")
        schema.addField("x_fp", "F", doc="x-coordinate of spot centroid in "
                        "focal plane coordinates (mm)")
        schema.addField("y_fp", "F", doc="y-coordinate of spot centroid in "
                        "focal plane coordinates (mm)")
        schema.addField("radius", "F", doc="width of spot in pixels")
        schema.addField("signal", "F", doc="Spot signal")
        schema.addField("bkg_mean", "F", doc="mean background")
        schema.addField("bkg_std", "F", doc="std of background")
        table = afw_table.SourceTable.make(schema)
        catalog = afw_table.SourceCatalog(table)
        for handle in exposure_handles:
            det = camera[handle.dataId['detector']]
            pix_to_fp = det.getTransform(cameraGeom.PIXELS,
                                         cameraGeom.FOCAL_PLANE)
            det_name = det.getName()
            dataId = handle.dataId
            exp_id, instrument = dataId['exposure'], dataId["instrument"]
            datasetType = handle.ref.datasetType.name
            collections = handle.ref.run
            idata = ImageData(repo=repo, collections=collections, instrument=instrument, dataId=dataId, datasetType=datasetType)
            if not forced:
                sp = Spot()
                image = idata.get_image()
                spot = sp.find_and_get_best_spots(image)
                ap = AperturePhotometry(image=idata, spot=spot)
                signal = ap.do_forced_aperture_photometry(centroid=ap.Spot.centroid, radius=aperture)
            elif forced :
                sp = Spot(x=ref_spots[handle.dataId['detector']]["x"], y=ref_spots[handle.dataId['detector']]["y"], radius=aperture)
                ap = AperturePhotometry(image=idata, spot=sp)
                signal = ap.do_forced_aperture_photometry(centroid=ap.Spot.centroid, radius=aperture)
            x, y = ap.Spot.x, ap.Spot.y
            radius = ap.Spot.radius
            fpx, fpy = pix_to_fp.getMapping().applyForward(np.vstack((x, y)))
            bkg_mean, bkg_std = ap.background_mean, ap.background_std
            if type(fpx)!=list:
                fpx, fpy=[fpx],[fpy]
            if type(x)!=list:
                x, y=[x],[y]
            if type(radius)!=list:
                radius=[radius]
            if type(signal)!=list:
                signal=[signal]
            for i in range(len(x)):
                record = catalog.addNew()
                record.set('detector', det_name)
                record.set('exposure', exp_id)
                record.set('x', x[i])
                record.set('y', y[i])
                record.set('x_fp', fpx[i])
                record.set('y_fp', fpy[i])
                record.set('radius', radius[i])
                record.set('signal', signal[i])
                record.set('bkg_mean', bkg_mean)
                record.set('bkg_std', bkg_std)
        return pipeBase.Struct(spot_catalog=catalog)
