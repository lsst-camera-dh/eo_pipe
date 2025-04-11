from lsst.afw import cameraGeom
import sys
import numpy as np
import lsst.afw.table as afw_table
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

__all__ = ["CBPSpotMeasurementTask"]


class CBPSpotMeasurementTaskConnections(pipeBase.PipelineTaskConnections,
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

    reference_spot_catalog = cT.PrerequisiteInput(
        name="cbp_spot_detection",
        doc="Catalog of cbp spot measurements.",
        storageClass="AstropyTable",
        dimensions=("instrument", "detector"))

    cbp_spot_detection = cT.Output(
        name="spotMeasurement",
        doc="Catalog of cbp spot measurements.",
        storageClass="SourceCatalog",
        dimensions=("instrument", "detector"))
    
    def __init__(self, *, config=None):
        super().__init__(config=config)

        if config.doForcedPhotometry is not True:
            self.prerequisiteInputs.remove("reference_spot_catalog")

class CBPSpotMeasurementTaskConfig(pipeBase.PipelineTaskConfig,
                                pipelineConnections=CBPSpotMeasurementTaskConnections):
    """Configuration for SpotMeasurementTask."""
    aperture = pexConfig.Field(doc="Aperture radius for spot photometry in pixels.",
                              dtype=int, default=200)  # 200 default, min=1 to avoid zero aperture

    #minarea = pexConfig.Field(doc="Minimum pixel area for spots.",
    #                          dtype=float, default=20000.0) # 30000 default
    #maxarea = pexConfig.Field(doc="Maximum pixel area for spots.",
    #""                          dtype=float, default=4000000.0)
    #force_circle = pexConfig.Field(doc="Force spot shape to be a circle",
    #                               dtype=bool, default=True)
    doForcedPhotometry = pexConfig.Field(
        dtype=bool,
        doc="Use forced photometry if True. If False, find new spots.",
        default=False,
    )

class CBPSpotMeasurementTask(pipeBase.PipelineTask):
    """Spot measurement task for CCOB narrow bean data."""
    ConfigClass = CBPSpotMeasurementTaskConfig
    _DefaultName = "cbpspotMeasurementTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, exposure_handles, camera,  reference_spot_catalog=None):
        sys.path.append("/sdf/home/a/amouroux/DATA/eo_pipe/python/lsst/eo/pipe")
        from photometry import AperturePhotometry, ImageData, Spot

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
        if self.config.doForcedPhotometry is True:
            print("Running forced photometry...")
            if reference_spot_catalog is None:
                raise RuntimeError("Forced photometry requested but no reference spot catalog provided.")
            ref_spots = {handle.dataId['detector']: reference_spot_catalog[handle.dataId['detector']] for handle in exposure_handles}
            if not ref_spots:
                raise RuntimeError("No reference spots found in the catalog.")
        for handle in exposure_handles:
            det = camera[handle.dataId['detector']]
            pix_to_fp = det.getTransform(cameraGeom.PIXELS,
                                         cameraGeom.FOCAL_PLANE)
            det_name = det.getName()
            #dataId = handle.dataId
            #exp_id, instrument = dataId['exposure'], dataId["instrument"]
            #datasetType = handle.ref.datasetType.name
            #collections = handle.ref.run
            exposure = handle.get()
            idata = ImageData(exposure_handle=handle)
            if self.config.doForcedPhotometry is False:
                sp = Spot()
                image = idata.get_image_from_handle()
                spot = sp.find_and_get_best_spot(idata.img.getImage())
                ap = AperturePhotometry(image=idata, spot=sp)
                signal = ap.do_forced_aperture_photometry(centroid=sp.centroid, radius=self.config.aperture)
            elif self.config.doForcedPhotometry is True :
                sp = Spot(x=ref_spots[handle.dataId['detector']]["x"], y=ref_spots[handle.dataId['detector']]["y"], radius=self.config.aperture)
                ap = AperturePhotometry(image=idata, spot=sp)
                signal = ap.do_forced_aperture_photometry(centroid=ap.Spot.centroid, radius=self.config.aperture)
            x, y = ap.Spot.x, ap.Spot.y
            radius = ap.Spot.radius
            if int(x)!=0:
                fpx, fpy = pix_to_fp.getMapping().applyForward(np.vstack((x, y)))
            else : 
                fpx, fpy = 0.0, 0.0
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
                record.set("exposure", handle.dataId['exposure'])
                record.set('detector', det_name)
                record.set('x', x[i])
                record.set('y', y[i])
                record.set('x_fp', fpx[i])
                record.set('y_fp', fpy[i])
                record.set('radius', radius[i])
                record.set('signal', signal[i])
                record.set('bkg_mean', bkg_mean)
                record.set('bkg_std', bkg_std)
        return pipeBase.Struct(cbp_spot_detection=catalog)
