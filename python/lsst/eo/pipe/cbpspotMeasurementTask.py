from lsst.afw import cameraGeom
import sys
import numpy as np
import lsst.afw.table as afw_table
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from astropy.table import Table, Column
import astropy.units as u

__all__ = ["CBPSpotMeasurementTask"]

class CBPSpotMeasurementTaskConnections(pipeBase.PipelineTaskConnections,
                                     dimensions=("instrument", "detector")):

    exposure_handles_input = cT.Input(
        name="post_isr_image",  
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

    reference_spot_catalog_input = cT.Input(
        name="cbp_spot_unforced",
        doc="Catalog of cbp spot measurements.",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "detector"))

    reference_spot_catalog_output = cT.Output(
        name="cbp_spot_unforced",
        doc="Catalog of cbp reference spots.",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "detector"))

    forced_spot_catalog_output = cT.Output(
        name="cbp_spot_forced",
        doc="Catalog of cbp spot measurements.",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "detector"))
    
    def __init__(self, *, config=None):
        super().__init__(config=config)
        if not config.doForcedPhotometry:
            del self.reference_spot_catalog_input
            del self.forced_spot_catalog_output
        else:
            del self.reference_spot_catalog_output

class CBPSpotMeasurementTaskConfig(pipeBase.PipelineTaskConfig,
                                pipelineConnections=CBPSpotMeasurementTaskConnections):
    """Configuration for SpotMeasurementTask."""
    aperture = pexConfig.Field(doc="Aperture radius for spot photometry in pixels.",
                              dtype=int, default=200)  # 200 default, min=1 to avoid zero aperture

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

    def run(self, exposure_handles, camera, reference_spot_catalog=None):
        sys.path.append("/sdf/home/a/amouroux/DATA/eo_pipe/python/lsst/eo/pipe")
        from photometry import AperturePhotometry, ImageData, Spot
        table = Table()
        if self.config.doForcedPhotometry is True:
            print("Running forced photometry...")
            if reference_spot_catalog is None:
                raise RuntimeError("Forced photometry requested but no reference spot catalog provided.")
            
        for handle in exposure_handles:
            det = camera[handle.dataId['detector']]
            pix_to_fp = det.getTransform(cameraGeom.PIXELS,
                                         cameraGeom.FOCAL_PLANE)
            det_name = det.getName()
            exposure = handle.get()
            print(f"Processing exposure {exposure}",
                 f"on detector {det_name} ; exposure_handle : {handle}")
            idata = ImageData(exposure_handle=handle)
            if self.config.doForcedPhotometry is False:
                sp = Spot(mask_size=100)
                print(f"idata exposure_handle:{idata.exposure_handle}")
                image = idata.get_image_from_handle()
                print(idata, image)
                spot = sp.find_and_get_best_spot(idata.img.getImage())
                print(f"spot : {spot}")
                print(idata.img.getImage(), idata.img, idata.img.getImage().getArray())
                ap = AperturePhotometry(image=idata, spot=sp)
                signal = ap.do_forced_aperture_photometry(centroid=sp.centroid, radius=self.config.aperture)
            elif self.config.doForcedPhotometry is True :
                ref_x, ref_y = np.median(reference_spot_catalog["x"]), np.median(reference_spot_catalog["y"])
                sp = Spot(x=ref_x, y=ref_y, radius=self.config.aperture)
                ap = AperturePhotometry(image=idata, spot=sp)
                signal = ap.do_forced_aperture_photometry(centroid=ap.Spot.centroid, radius=self.config.aperture)
            x, y = ap.Spot.x, ap.Spot.y
            radius = ap.Spot.radius
            exposure_time = idata.shuttime 
            obs_annot = idata.obsannot
            try :
                wavelength = int(obs_annot)
            except :
                wavelength = 0
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
            table["exposure"] = Column([handle.dataId['exposure']] * len(x), unit=u.dimensionless_unscaled, description="Exposure ID")
            table["detector"] = Column([det_name] * len(x), unit=u.dimensionless_unscaled, description="Detector name")
            table["x"] = Column(x, unit=u.pixel, description="X-coordinate of spot centroid on CCD")
            table["y"] = Column(y, unit=u.pixel, description="Y-coordinate of spot centroid on CCD")
            table["x_fp"] = Column(fpx, unit=u.mm, description="X-coordinate in focal plane coordinates")
            table["y_fp"] = Column(fpy, unit=u.mm, description="Y-coordinate in focal plane coordinates")
            table["signal"] = Column(signal, unit=u.adu, description="Integrated signal of the spot")
            table["radius"] = Column(radius, unit=u.pixel, description="Radius of the spot or aperture if forced")
            table["bkg_mean"] = Column(bkg_mean, unit=u.adu, description="Mean background level")
            table["bkg_std"] = Column(bkg_std, unit=u.adu, description="Standard deviation of background level")
            table["exposure_time"] = Column([exposure_time] * len(x), unit=u.s, description="Exposure time")
            table["wavelength"] = Column([wavelength] * len(x), unit=u.nm, description="Wavelength of observation")
        return pipeBase.Struct(cbp_spot_detection=table)