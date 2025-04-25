import numpy as np
from astropy.io import fits
from astropy.table import Table
from lsst.daf.butler import Butler
import lsst.afw.detection as afwDetect
from lsst.afw.geom import SpanSet
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry, ApertureStats, CircularAnnulus
from photutils.background import Background2D
import sys, os

sys.path.append('/sdf/group/rubin/user/amouroux/comissioning/cbp_analysis/python/lsst/python/exposure_time_calculator')
from rubin_calib_etc import RubinCalibETC

class Spot:
    def __init__(self, x=None, y=None, radius=None, mask_size = 150): 
        self.magnification_factor = 10.31 / 0.635
        self.mask_size = mask_size #mask size in um
        self.centroid = (x,y)
        self.x = x
        self.y = y
        self.radius = radius
        if self.x is not None and self.y is not None:
            self.centroid = (self.x, self.y)

    def get_mask_size(self):
        self.mask_size_fp = self.mask_size * self.magnification_factor
        self.mask_size_fp_px = self.mask_size_fp / 10 # in px 1px = 10um
        self.mask_area_fp_px = np.pi * (self.mask_size_fp_px / 2) ** 2

    def find_spot(self, image, threshold_adu=100, minarea=20000):
        """
        Find spots in the image using a threshold and minimum area.

        Parameters:
        threshold_adu (int): Threshold in ADU.
        minarea (int): Minimum area of the spot.
        """
        threshold = afwDetect.Threshold(threshold_adu)
        self.found_spot = afwDetect.FootprintSet(image, threshold, npixMin=minarea).getFootprints()
        return self.found_spot

    def find_spot_iteratively(self, image, threshold_adu_min=100, threshold_adu_max=500, minarea=20000):
        """
        Find spots in the image using a threshold and minimum area.

        Parameters:
        threshold_adu (int): Threshold in ADU.
        minarea (int): Minimum area of the spot.
        """
        threshold_array = np.logspace(np.log10(threshold_adu_max), np.log10(threshold_adu_min),5)
        for threshold_adu in threshold_array:
            found_spot = self.find_spot(image, threshold_adu, minarea)
            if len(found_spot) > 0:
                break
        self.found_spot = found_spot
        return self.found_spot
    
    def get_spot_information(self, spot=None):
        """
        Get the centroid and radius of a spot.

        Parameters:
        spot: The spot to calculate the centroid and radius for.

        Returns:
        tuple: Centroid and radius of the spot.
        """
        if spot is None:
            spot = self.found_spot[0]
        self.centroid = spot.getCentroid()
        self.x = self.centroid.getX()
        self.y = self.centroid.getY()
        self.radius = np.sqrt(spot.getArea() / np.pi)

    def get_best_spot_radius(self, spots=None):
        """
        Get the best spot from a list of spots.

        Parameters:
        spots (list): List of spots to choose from.

        Returns:
        Spot: The best spot.
        """
        if spots is None:
            spots = self.found_spot
        spots.sort(key=lambda x: x.getArea(), reverse=True)
        self.get_mask_size()
        if len(spots) > 1:
            print("Multiple spots found, choosing the largest one")
        for spot in spots:
            radius = np.sqrt(spot.getArea() / np.pi)
            if radius > self.mask_size_fp_px*2:
                print("Spot found with radius more than 2 times larger than mask size. Passing to the next spot")
                continue
            elif radius < self.mask_size_fp_px/4:
                print("Spot found with radius less than half of mask size. Passing to the next spot")
                continue
            else:
                print("Best spot found")
                self.best_spot = spot
                break
        return self.best_spot

    def get_best_spot(self, image, spots=None):
        """
        Find the brightest spot from a list of spots.

        Parameters:
        image: The image to analyze.
        spots (list): List of spots to analyze. If None, spots will be found in the image.

        Returns:
        Spot: The brightest spot.
        """
        if spots is None:
            spots = self.find_spot(image)

        spot_counts = []
        for spot in spots:
            self.get_spot_information(spot)
            ap = AperturePhotometry(image=None, spot=spot)
            aperture = ap.generate_aperture(centroid=(spot.x, spot.y), radius=200)
            count = ap.do_aperture_photometry(image=image, aperture=aperture)
            spot_counts.append((spot, count))

        brightest_spot = max(spot_counts, key=lambda x: x[1])[0]
        return brightest_spot
    
    def find_and_get_best_spot(self, image):
        self.get_mask_size()
        self.find_spot_iteratively(image, threshold_adu_min=100, threshold_adu_max=1e4, minarea=int(self.mask_area_fp_px*.9))
        if len(self.found_spot) == 0:
            print("No spot found")
            ss = SpanSet.fromShape(0, offset=(0,0))
            self.best_spot = afwDetect.Footprint(ss)
        elif len(self.found_spot) > 1:
            self.get_best_spot(image, spots=self.found_spot)
        else:
            self.best_spot = self.found_spot[0]
        self.get_spot_information(self.best_spot)
        return self.best_spot

class ImageData:
    def __init__(self, exposure_handle=None, repo=None):
        self.repo = repo
        self.exposure_handle = exposure_handle
        if exposure_handle is not None and repo is not None:
            self.dataId = exposure_handle.dataId
            self.instrument= self.dataId["instrument"]
            self.datasetType = exposure_handle.ref.datasetType.name
            self.collections = exposure_handle.ref.run
        self.image = None
        return None
    
    def get_datasets(self):
        self.butler = Butler(self.repo)
        self.registry = self.butler.registry
        self.datasets = list(self.registry.queryDatasets(datasetType=self.datasetType, collections=self.collections, instrument=self.instrument, dataId=self.dataId))
        if len(self.datasets) == 0:
            print("No dataset found")
        elif len(self.datasets) > 1:
            print("Multiple datasets found")
        elif len(self.datasets) == 1:
            print("One dataset found")
        return self.datasets
    
    def get_image(self, dataset=None):
        if dataset is None:
            dataset = self.get_datasets()[0]
        self.img = self.butler.get(dataset)
        self.metadata = dict(self.img.getMetadata())
        self.image = self.img.getImage().getArray()
        return self.image

    def get_image_from_handle(self, exposure_handle=None):
        """
        Get the image from a given exposure object.

        Parameters:
        exposure: The exposure object to get the image from.

        Returns:
        np.ndarray: The image array.
        """
        if exposure_handle is None:
            exposure_handle = self.exposure_handle
        self.img = self.exposure_handle.get()
        self.metadata = dict(self.img.getMetadata())
        self.image = self.img.getImage().getArray()
        return self.image
    
class AperturePhotometry:
    """
    A class to perform aperture photometry on CBP images.
    """
    def __init__(self, image=None, spot=None):
        self.ImageData = image #ImageData object 
        self.Spot = spot #Spot object
        if self.ImageData is not None:
            if self.ImageData.image is None and self.ImageData.exposure_handle is None:
                self.ImageData.get_image()
            elif self.ImageData.image is None and self.ImageData.exposure_handle is not None:
                self.ImageData.get_image_from_handle()
            self.image = self.ImageData.image
        self.ImageData.shuttime = self.ImageData.metadata["SHUTTIME"]
        self.ImageData.obsannot = self.ImageData.metadata["OBSANNOT"]
        self.background = None

    def get_2d_background_threshold(self, threshold = None):
        """
        Calculate the background of the image.
        """
        if threshold is None:
            threshold = np.mean(self.image) + (3 * np.std(self.image))
        mask = np.ma.masked_where(self.image > threshold, self.image)
        bkg = Background2D(self.image, (int(len(self.image)/10), int(len(self.image[0])/10)), mask=mask.mask, exclude_percentile=50.0)
        self.background = bkg.background
        return self.background
    
    def get_2d_background_aperture(self, aperture = None):
        """
        Calculate the background of the image.
        """
        if aperture is None:
            aperture = self.aperture
        mask = np.zeros((self.image.shape[0], self.image.shape[1]), dtype=bool)
        mask |= aperture.to_mask(method='center').to_image((self.image.shape[0], self.image.shape[1])).astype(bool)
        bkg = Background2D(self.image, (int(len(self.image)/10), int(len(self.image[0])/10)), mask=mask, exclude_percentile=50.0)
        self.background = bkg.background
        return self.background

    def get_mean_background(self, threshold = None):
        """
        Calculate the background of the image.
        """
        if threshold is None:
            threshold = np.mean(self.image) + (3 * np.std(self.image))
        mask = np.zeros((self.image.shape[0], self.image.shape[1]), dtype=bool)
        mask |= (self.image > threshold)
        self.background = np.mean(self.image[mask==False])
        return self.background
    
    def get_dark_background(self): #To change with the task so that this is an input parameter ex : doBackground substraction | background_substraction_method
        exp_table_path = f"/sdf/group/rubin/user/amouroux/comissioning/cbp_analysis/notebooks/comcam_analysis/exposures/{self.ImageData.img.metadata['FILTBAND']}_{self.ImageData.img.metadata['DAYOBS']}.fits"
        exp_table = Table.read(exp_table_path)
        dark_exp = exp_table[exp_table["exposure"] == self.ImageData.dataId["exposure"]]["dark_exposure"][0]
        dataId_dark = {"exposure": dark_exp, "detector": self.ImageData.dataId["detector"]}
        image_data_copy = self.ImageData
        image_data_copy.dataId = dataId_dark
        dark_image = image_data_copy.get_image()
        self.background = dark_image
        return self.background
    
    def get_annulus_background(self, annulus=None):
        """
        Calculate the background of the image.
        """
        if annulus is None:
            annulus = CircularAnnulus(self.Spot.centroid, self.Spot.radius + 100, self.Spot.radius + 300)
        self.annulus_stats = ApertureStats(self.image, annulus)
        self.background = self.annulus_stats.mean
        return self.background
    
    def get_substracted_background_image(self, background=None):
        """
        Subtract the background from the image.
        """
        if background is None:
            background = self.background
            if self.background is None :
                print("No background found")
                return None
        self.subtracted_background_image = self.image - background 
        return self.subtracted_background_image
    
    def generate_aperture(self, centroid=None, radius=None):
        """
        Generate an aperture.
        """
        if centroid is None:
            centroid = (self.Spot.x, self.Spot.y)
        if radius is None:
            radius = self.Spot.radius
        self.aperture = CircularAperture(centroid, r=radius)
        return self.aperture
    
    def do_aperture_photometry(self, image = None, aperture=None):
        """
        Perform aperture photometry on the image.
        """
        if aperture is None:
            aperture = self.aperture # aperture object
        if image is None:
            image = self.subtracted_background_image
        self.adu_count = aperture.do_photometry(image)[0][0]
        return self.adu_count
    
    def do_forced_aperture_photometry(self, centroid=None, radius=None):
        """
        Perform aperture photometry on the image.
        """
        if centroid is None:
            centroid = (self.Spot.x, self.Spot.y)
        if radius is None:
            radius = self.Spot.radius
        background_aperture = CircularAperture(centroid, r=600)
        background = self.get_2d_background_aperture(aperture=background_aperture)
        substracted_background_image = self.get_substracted_background_image(background=background)
        aperture = self.generate_aperture(centroid=centroid, radius=radius)
        self.background_stats = ApertureStats(background, aperture)
        self.background_mean, self.background_std = self.background_stats.mean, self.background_stats.std
        adu_count = self.do_aperture_photometry(image = substracted_background_image, aperture = aperture)
        return adu_count

    ######### Useless for now ##########    
    def get_detector_efficiency(self, wavelength, exptime_config = "/sdf/data/rubin/user/amouroux/comissioning/cbp_codes/notebooks/imsim/cbp_calib.yaml"):
        try:
            ETC = RubinCalibETC(exptime_config)
            ETC.get_detector_efficiency()
            d_eff = ETC.detector_efficiency
            w_index = np.where(ETC.rubin_wavelength == float(wavelength * 1e9))[0][0]
            self.detector_efficiency = d_eff[w_index]
        except Exception as e:
            print(f"Error calculating detector efficiency: {e}")
            self.detector_efficiency = None
    
    def adu_to_electron(self, gain=1):
        """
        Convert ADU counts to electron counts.

        Parameters:
        gain (float): Gain factor to convert ADU to electrons.

        Returns:
        np.ndarray: Electron counts.
        """
        self.el_count = np.array(self.adu_count) * gain
        return self.el_count

    def electron_to_photon(self):
        """
        Convert electron counts to photon counts.

        Returns:
        np.ndarray: Photon counts.
        """
        if self.detector_efficiency is None:
            raise ValueError("Detector efficiency is not set. Please run get_detector_efficiency first.")
        
        h = 6.62607015e-34  # Planck constant in J*s
        c = 3e8  # Speed of light in m/s
        self.n_photons = self.el_count / self.detector_efficiency
        return self.n_photons

    def get_flux(self, count, total_exptime=13):
        """
        Calculate the flux from the count and exposure time.

        Parameters:
        count (float): Count value.
        total_exptime (float): Total exposure time in seconds.

        Returns:
        float: Flux value.
        """
        area = self.found_spot.getArea() * (1e-5)**2  # Area in m^2
        self.flux = count / (area * total_exptime)
        return self.flux

    def get_physical_photon_flux(self, total_exptime=13):
        """
        Calculate the physical photon flux.

        Parameters:
        total_exptime (float): Total exposure time in seconds.

        Returns:
        float: Physical photon flux in W/m^2.
        """
        if self.wavelength is None:
            raise ValueError("Wavelength is not set. Please set the wavelength first.")
        
        h = 6.62607015e-34  # Planck constant in J*s
        c = 3e8  # Speed of light in m/s
        E_photon = h * c / self.wavelength  # Energy per photon in J
        E_tot = self.n_photons * E_photon  # Total energy in J
        area = self.found_spot.getArea() * (1e-5)**2  # Area in m^2
        self.physical_photon_flux = E_tot / (area * total_exptime)  # Flux in W/m^2
        return self.physical_photon_flux