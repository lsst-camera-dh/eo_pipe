import numpy as np
from astropy.table import Table
from lsst.daf.butler import Butler
import lsst.afw.detection as afwDetect
from lsst.afw.geom import SpanSet
from photutils.aperture import CircularAperture, ApertureStats, CircularAnnulus
from photutils.background import Background2D#, BkgZoomInterpolator

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

    def _check_image_quality(self, image, min_signal_ratio=10.0):
        """
        Check if the image has sufficient signal to contain a spot.
        
        Parameters:
        min_signal_ratio: Minimum ratio between max pixel and median
        """
        image_array = image.array
        median_val = np.median(image_array)
        max_val = np.percentile(image_array, 99.99)
        std_val = np.std(image_array)
        
        signal_ratio = (max_val - median_val) / std_val
        print(signal_ratio)
        if signal_ratio < min_signal_ratio:
            print(f"Low signal detected (ratio: {signal_ratio:.2f}). Skipping spot detection.")
            return False
        return True

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

    def find_spot_iteratively(self, image, threshold_array, minarea=20000):
        """
        Find spots in the image using a threshold and minimum area.

        Parameters:
        threshold_adu (int): Threshold in ADU.
        minarea (int): Minimum area of the spot.
        """
        #threshold_array = np.logspace(np.log10(threshold_adu_max), np.log10(threshold_adu_min),5)
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
            self.get_spot_information(spot=spot)
            ap = AperturePhotometry(image=None, spot=None)
            aperture = ap.generate_aperture(centroid=(self.x, self.y), radius=200)
            count = ap.do_aperture_photometry(image=image.array, aperture=aperture)
            spot_counts.append((spot, count))

        brightest_spot = max(spot_counts, key=lambda x: x[1])[0]
        return brightest_spot
    
    def find_and_get_best_spot(self, image):
        self.get_mask_size()
        threshold_array = np.array([ np.median(image.array)+(25 * np.std(image.array)),
                            np.median(image.array)+(5 * np.std(image.array)),
                            np.percentile(image.array, 95),
                            np.percentile(image.array, 75),
                            np.percentile(image.array, 60),
                            np.median(image.array)
                          ])
        if not self._check_image_quality(image):
            self.found_spot = []
        else:
            self.find_spot_iteratively(image, threshold_array=threshold_array, minarea=int(self.mask_area_fp_px*.9))
        if len(self.found_spot) == 0:
            print("No spot found")
            ss = SpanSet.fromShape(0, offset=(0,0))
            self.best_spot = afwDetect.Footprint(ss)
        elif len(self.found_spot) > 1:
            print(f"Found {len(self.found_spot)} spots.")
            self.best_spot = self.get_best_spot(image, spots=self.found_spot)
        else:
            print("Found one spot.")
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
        """   
        interpolator = BkgZoomInterpolator(
            order=3,          # cubic interpolation (0=nearest, 1=linear, 3=cubic, etc.)
            mode='reflect',   # how to handle edges ('reflect', 'constant', 'nearest', etc.)
            cval=0.0,         # used if mode='constant'
            clip=True,        # clip interpolated values to min/max of original boxes
            grid_mode=True    # align interpolation grid with input image pixels
            )
        """
        mask = np.zeros((self.image.shape[0], self.image.shape[1]), dtype=bool)
        mask |= aperture.to_mask(method='center').to_image((self.image.shape[0], self.image.shape[1])).astype(bool)
        bkg = Background2D(self.image, (int(len(self.image)/10), int(len(self.image[0])/10)), mask=mask, exclude_percentile=50.0)
        #bkg = Background2D(self.image, (int(len(self.image)/40), int(len(self.image)/40)), mask=mask, exclude_percentile=20.0,interpolator = interpolator)
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
        Add do background_substraction
        """
        if centroid is None:
            centroid = (self.Spot.x, self.Spot.y)
        if radius is None:
            radius = self.Spot.radius
        background_aperture = CircularAperture(centroid, r=600)
        background = self.get_2d_background_aperture(aperture=background_aperture)
        substracted_background_image = self.get_substracted_background_image(background=background)
        #substracted_background_image = self.image
        aperture = self.generate_aperture(centroid=centroid, radius=radius)
        self.background_stats = ApertureStats(background, aperture)
        self.background_mean, self.background_std = self.background_stats.mean, self.background_stats.std
        adu_count = self.do_aperture_photometry(image = substracted_background_image, aperture = aperture)
        return adu_count