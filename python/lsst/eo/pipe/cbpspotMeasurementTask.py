from lsst.afw import cameraGeom
import os, sys
import numpy as np
import lsst.afw.table as afw_table
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from astropy.table import Table, Column
import astropy.units as u
import cProfile
sys.path.append(os.path.dirname(__file__))
from photometry import AperturePhotometry, ImageData, Spot

__all__ = ["CBPSpotMeasurementTask"]

class CBPSpotMeasurementTaskConnections(pipeBase.PipelineTaskConnections,
                                     dimensions=("instrument", "detector")):

    exposure_handles = cT.Input(
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

    doMultipleDetections = pexConfig.Field(
        dtype=bool,
        doc="Allows to detect several spots on the same exposure. For reference spot mode only.",
        default=False,
    )

    doMultipleApertures = pexConfig.Field(
        dtype=bool,
        doc="Allows to use multiple apertures for forced photometry.",
        default=False,
    )

class CBPSpotMeasurementTask(pipeBase.PipelineTask):
    """Spot measurement task for CCOB narrow bean data."""
    ConfigClass = CBPSpotMeasurementTaskConfig
    _DefaultName = "cbpspotMeasurementTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def forced_photometry_multiple_apertures(self, exposure_handles, camera, reference_spot_catalog_input=None):
        """Perform forced photometry on the reference spots."""
        data_dict = {
            "exposure": [],
            "detector": [],
            "x": [],
            "y": [],
            "x_fp": [],
            "y_fp": [],
            "signal": [],
            "radius": [],
            "bkg_mean": [],
            "bkg_std": [],
            "exposure_time": [],
            "wavelength": [],
        }

        if reference_spot_catalog_input is None:
                raise RuntimeError("Forced photometry requested but no reference spot catalog provided.")
        else :
            if len(reference_spot_catalog_input) == 0:
                raise RuntimeError("Reference spot catalog is empty.")
            else :
                print(f"Reference spot catalog provided : \n{reference_spot_catalog_input}\n")
                print("Getting brightest spot from reference catalog...\n")
                reference_spot_catalog_input.sort("signal", reverse=True)
                reference_spots = reference_spot_catalog_input[:3]
                print(f"Reference spots : \n{reference_spots}\n")

        det = camera[exposure_handles[0].dataId['detector']]
        pix_to_fp = det.getTransform(cameraGeom.PIXELS,
                                         cameraGeom.FOCAL_PLANE)
        det_name = det.getName()
        ref_x, ref_y = np.median(reference_spots["x"]), np.median(reference_spots["y"])
        sp = Spot(x=ref_x, y=ref_y, radius=self.config.aperture)
        fpx, fpy = pix_to_fp.getMapping().applyForward(np.vstack((ref_x, ref_y)))
        for handle in exposure_handles:
            exposure = handle.dataId["exposure"]
            print(f"Processing exposure {exposure}",
                f"on detector {det_name}")
            idata = ImageData(exposure_handle=handle)
            ap = AperturePhotometry(image=idata, spot=sp)
            aperture_radii = np.array([10,20,30,40,50,75,100,150,200,250,300,350,400,450,500])
            signal_list = []
            for aperture_radius in aperture_radii:
                signal_list.append(ap.do_forced_aperture_photometry(centroid=ap.Spot.centroid, radius=aperture_radius))
            exposure_time = ap.ImageData.shuttime 
            obs_annot = ap.ImageData.obsannot
            print("exposure_time : ", exposure_time, "obs_annot : ", obs_annot)
            try:
                wavelength = int(obs_annot)
            except :
                wavelength = 0
            for i, signal in enumerate(signal_list):
                data_dict["exposure"].append(exposure)
                data_dict["detector"].append(det_name)
                data_dict["exposure_time"].append(exposure_time)
                data_dict["wavelength"].append(wavelength)
                data_dict["x"].append(ref_x)
                data_dict["y"].append(ref_y)
                data_dict["x_fp"].append(fpx[0])
                data_dict["y_fp"].append(fpy[0])
                data_dict["signal"].append(signal)
                data_dict["radius"].append(aperture_radii[i])
                data_dict["bkg_mean"].append(ap.background_mean)
                data_dict["bkg_std"].append(ap.background_std)
        return data_dict
  
    def forced_photometry(self, exposure_handles, camera, reference_spot_catalog_input=None):
        """Perform forced photometry on the reference spots."""
        data_dict = {
            "exposure": [],
            "detector": [],
            "x": [],
            "y": [],
            "x_fp": [],
            "y_fp": [],
            "signal": [],
            "radius": [],
            "bkg_mean": [],
            "bkg_std": [],
            "exposure_time": [],
            "wavelength": [],
        }

        if reference_spot_catalog_input is None:
                raise RuntimeError("Forced photometry requested but no reference spot catalog provided.")
        else :
            if len(reference_spot_catalog_input) == 0:
                raise RuntimeError("Reference spot catalog is empty.")
            else :
                print(f"Reference spot catalog provided : \n{reference_spot_catalog_input}\n")
                print("Getting brightest spot from reference catalog...\n")
                reference_spot_catalog_input.sort("signal", reverse=True)
                reference_spots = reference_spot_catalog_input[:3]
                print(f"Reference spots : \n{reference_spots}\n")

        det = camera[exposure_handles[0].dataId['detector']]
        pix_to_fp = det.getTransform(cameraGeom.PIXELS,
                                         cameraGeom.FOCAL_PLANE)
        det_name = det.getName()
        ref_x, ref_y = np.median(reference_spots["x"]), np.median(reference_spots["y"])
        sp = Spot(x=ref_x, y=ref_y, radius=self.config.aperture)
        fpx, fpy = pix_to_fp.getMapping().applyForward(np.vstack((ref_x, ref_y)))
        for handle in exposure_handles:
            exposure = handle.dataId["exposure"]
            print(f"Processing exposure {exposure}",
                f"on detector {det_name}")
            idata = ImageData(exposure_handle=handle)
            ap = AperturePhotometry(image=idata, spot=sp)
            ####TEMP : growth of curve
            signal = ap.do_forced_aperture_photometry(centroid=ap.Spot.centroid, radius=self.config.aperture)
            exposure_time = ap.ImageData.shuttime 
            obs_annot = ap.ImageData.obsannot
            print("exposure_time : ", exposure_time, "obs_annot : ", obs_annot)
            try :
                wavelength = int(obs_annot)
            except :
                wavelength = 0
            data_dict["exposure"].append(exposure)
            data_dict["detector"].append(det_name)
            data_dict["exposure_time"].append(exposure_time)
            data_dict["wavelength"].append(wavelength)
            data_dict["x"].append(ref_x)
            data_dict["y"].append(ref_y)
            data_dict["x_fp"].append(fpx[0])
            data_dict["y_fp"].append(fpy[0])
            data_dict["signal"].append(signal)
            data_dict["radius"].append(self.config.aperture)
            data_dict["bkg_mean"].append(ap.background_mean)
            data_dict["bkg_std"].append(ap.background_std)
        return data_dict
        
    def unforced_photometry_single(self, exposure_handles, camera):
        """Perform unforced photometry on a single spot."""
        data_dict = {
            "exposure": [],
            "detector": [],
            "x": [],
            "y": [],
            "x_fp": [],
            "y_fp": [],
            "signal": [],
            "radius": [],
            "bkg_mean": [],
            "bkg_std": [],
            "exposure_time": [],
            "wavelength": [],
        }
        det = camera[exposure_handles[0].dataId['detector']]
        pix_to_fp = det.getTransform(cameraGeom.PIXELS,
                                         cameraGeom.FOCAL_PLANE)
        det_name = det.getName()
        for handle in exposure_handles:
            exposure = handle.dataId["exposure"]
            print(f"Processing exposure {exposure}",
                f"on detector {det_name}")
            sp = Spot(mask_size=100)
            idata = ImageData(exposure_handle=handle)
            image = idata.get_image_from_handle()
            print("Finding best spot...")
            spot = sp.find_and_get_best_spot(idata.img.getImage())
            print(f"Spot : {spot}")
            centroid = spot.getCentroid()
            x_sp, y_sp = centroid.getX(), centroid.getY()
            radius = np.sqrt(spot.getArea() / np.pi)
            sp = Spot(x=x_sp, y=y_sp, radius=radius)
            fpx, fpy = pix_to_fp.getMapping().applyForward(np.vstack((x_sp, y_sp)))
            ap = AperturePhotometry(image=idata, spot=sp)
            signal = ap.do_forced_aperture_photometry(centroid=sp.centroid, radius=self.config.aperture)
            exposure_time = ap.ImageData.shuttime 
            obs_annot = ap.ImageData.obsannot
            print("exposure_time : ", exposure_time, "obs_annot : ", obs_annot)
            try :
                wavelength = int(obs_annot)
            except :
                wavelength = 0
            data_dict["exposure"].append(exposure)
            data_dict["detector"].append(det_name)
            data_dict["exposure_time"].append(exposure_time)
            data_dict["wavelength"].append(wavelength)
            data_dict["x"].append(x_sp)
            data_dict["y"].append(y_sp)
            data_dict["x_fp"].append(fpx[0])
            data_dict["y_fp"].append(fpy[0])
            data_dict["signal"].append(signal)
            data_dict["radius"].append(radius)
            data_dict["bkg_mean"].append(ap.background_mean)
            data_dict["bkg_std"].append(ap.background_std)
        return data_dict

    def unforced_photometry_multiple(self, exposure_handles, camera):
        """Perform unforced photometry on multiple spots."""
        data_dict = {
            "exposure": [],
            "detector": [],
            "x": [],
            "y": [],
            "x_fp": [],
            "y_fp": [],
            "signal": [],
            "radius": [],
            "bkg_mean": [],
            "bkg_std": [],
            "exposure_time": [],
            "wavelength": [],
        }
        det = camera[exposure_handles[0].dataId['detector']]
        pix_to_fp = det.getTransform(cameraGeom.PIXELS,
                                         cameraGeom.FOCAL_PLANE)
        det_name = det.getName()
        for handle in exposure_handles:
            exposure = handle.dataId["exposure"]
            print(f"Processing exposure {exposure}",
                f"on detector {det_name}")
            idata = ImageData(exposure_handle=handle)
            image = idata.get_image_from_handle()
            sp = Spot(mask_size=10)
            threshold_array = np.array([ np.median(image)+(25 * np.std(image)),
                    np.median(image)+(5 * np.std(image))])
            spots = sp.find_spot_iteratively(idata.img.getImage(), threshold_array,minarea=200)
            if len(spots) > 25:
                print("Too many spots detected, keeping only the 25 largest ones...")
                spots = sorted(spots, key=lambda s: s.getArea(), reverse=True)[:25]
            for s in spots:
                centroid = s.getCentroid()
                x_sp, y_sp = centroid.getX(), centroid.getY()
                radius = np.sqrt(s.getArea() / np.pi)
                sp = Spot(x=x_sp, y=y_sp, radius=radius)
                fpx, fpy = pix_to_fp.getMapping().applyForward(np.vstack((x_sp, y_sp)))
                ap = AperturePhotometry(image=idata, spot=sp)
                aperture = ap.generate_aperture(centroid = centroid, radius=30)
                signal = ap.do_aperture_photometry(image=image, aperture=aperture)
                exposure_time = ap.ImageData.shuttime 
                obs_annot = ap.ImageData.obsannot
                print("exposure_time : ", exposure_time, "obs_annot : ", obs_annot)
                try :
                    wavelength = int(obs_annot)
                except :
                    wavelength = 0
                data_dict["exposure"].append(exposure)
                data_dict["detector"].append(det_name)
                data_dict["exposure_time"].append(exposure_time)
                data_dict["wavelength"].append(wavelength)
                data_dict["x"].append(x_sp)
                data_dict["y"].append(y_sp)
                data_dict["x_fp"].append(fpx[0])
                data_dict["y_fp"].append(fpy[0])
                data_dict["signal"].append(signal)
                data_dict["radius"].append(radius)
                data_dict["bkg_mean"].append(0.0)
                data_dict["bkg_std"].append(0.0)
        return data_dict
    

    def run(self, exposure_handles, camera, reference_spot_catalog_input=None):
        profiler = cProfile.Profile()
        profiler.enable() 
        print(f"Number of exposures to process: {len(exposure_handles)}")

        if not self.config.doForcedPhotometry :
            if not self.config.doMultipleDetections:
                print("Running unforced photometry for single spot...\n")
                data_dict = self.unforced_photometry_single(exposure_handles, camera)
            else: 
                print("Running unforced photometry for multiple spots...\n")
                data_dict = self.unforced_photometry_multiple(exposure_handles, camera)

        else : # do forced photometry change here for multiple apertures
            print("Running forced photometry...\n")
            if not self.config.doMultipleApertures:
                data_dict = self.forced_photometry(exposure_handles, camera, reference_spot_catalog_input)
            else:
                data_dict = self.forced_photometry_multiple_apertures(exposure_handles, camera, reference_spot_catalog_input)

        print("Creating table...")
        table = Table()
        unit_desc_dict = {
            "exposure": (u.dimensionless_unscaled, "Exposure ID"),
            "detector": (u.dimensionless_unscaled, "Detector name"),
            "exposure_time": (u.s, "Exposure time in seconds"),
            "wavelength": (u.nm, "Wavelength of observation in nanometers"),
            "x": (u.pixel, "X-coordinate of spot centroid on CCD"),
            "y": (u.pixel, "Y-coordinate of spot centroid on CCD"),
            "x_fp": (u.mm, "X-coordinate in focal plane coordinates"),
            "y_fp": (u.mm, "Y-coordinate in focal plane coordinates"),
            "signal": (u.adu, "Integrated signal of the spot"),
            "radius": (u.pixel, "Radius of the spot or aperture if forced"),
            "bkg_mean": (u.adu, "Mean background level"),
            "bkg_std": (u.adu, "Standard deviation of background level"),
        }
        for key, (unit, description) in unit_desc_dict.items():
            table[key] = Column(data_dict[key], unit=unit, description=description)

        print(f"Final table summary: {table.info}")
        print(f"First few rows: \n{table[:5]}")
        profiler.disable() 
        profiler.dump_stats("/sdf/home/a/amouroux/rubin-user/stack/w_2025_18/eo_pipe/new_task/cbpspotMeasurementTask_modifs.prof")
        print("Profiling data saved to cbpspotMeasurementTask.prof")
        if not self.config.doForcedPhotometry:
            return pipeBase.Struct(reference_spot_catalog_output=table)
        else:
            return pipeBase.Struct(forced_spot_catalog_output=table)