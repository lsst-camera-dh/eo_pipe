import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from lsst.daf.butler import Butler
import lsst.afw.detection as afwDetect
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry, ApertureStats
from photutils.background import Background2D
from matplotlib.patches import Circle
import sys
import yaml
sys.path.append('/sdf/group/rubin/user/amouroux/comissioning/cbp_analysis/python/lsst/python/exposure_time_calculator')
from rubin_calib_etc import RubinCalibETC

class AperturePhotometry():
    def __init__(self, config_file= "../../../../data/photometry_config.yaml"):
        with open(config_file) as f:
            config = yaml.safe_load(f)
            #print(config)
        self.repo_path = config["repo"]
        self.collections = config["collections"]
        self.instrument = config["instrument"]
        self.detector = config["detector"]
        self.datasetType = config["datasetType"]
        self.exposure=config["exposure"]
        self.spot_size = config["spot_size"]
        self.exptime_config = config["exptime_config"]
        self.wavelength = config["wavelength"]
        self.filter = config["filter"]
        self.total_exptime = config["exptime"]
        self.magnification_factor = 10.31 / 0.635
        return None
        
    def get_image(self):
        butler = Butler(self.repo_path,collections=self.collections)
        registry = butler.registry
        datasets = list(registry.queryDatasets(self.datasetType, instrument=self.instrument, detector=self.detector,
                                               exposure=self.exposure, collections=self.collections))
        #print(datasets)
        self.img = butler.get(datasets[0])
        self.image = self.img.image.array
    
    def get_spot_area(self):
        spot_size_fp = self.spot_size * self.magnification_factor
        self.spot_size_fp_px = spot_size_fp/1e-5
        self.spot_area = np.pi*(self.spot_size_fp_px**2)
        #print("Spot size on focal plane[px]=", self.spot_size_fp_px,"\n", "Spot area[px^2]=",self.spot_area)
    
    def find_spot(self, threshold_adu=100, minarea=20000):
        threshold = afwDetect.Threshold(threshold_adu)
        self.found_spot = afwDetect.FootprintSet(self.img.getImage(), threshold, npixMin=minarea).getFootprints()
    
    def get_image_mu_sig(self):
        self.mean, self.std = np.mean(self.image), np.std(self.image)
    
    def get_centroid(self, spot):
        centroid = spot.getCentroid()
        radius = np.sqrt(spot.getArea()/np.pi)
        #print("Spot centroid coordinate :", centroid, "\nSpot radius :", radius)
        return centroid, radius
        
    def get_background(self):
        mask = np.ma.masked_where(self.image>(self.mean+5*self.std), self.image)  
        bkg = Background2D(self.image, (int(len(self.image)/10), int(len(self.image[0])/10)), mask = mask.mask)
        self.background = bkg.background

    def get_substracted_background_image(self):
        self.get_background()
        self.substracted_background_image = self.image-self.background

    def do_aperture_photometry(self):
        self.get_image()
        self.get_spot_area()
        self.get_image_mu_sig()
        self.find_spot(threshold_adu=self.mean+5*self.std, minarea = int(.9*self.spot_area))
        n_spot = len(self.found_spot)
        if n_spot>1:
            self.centroids, self.radiuses = [], []
            self.x, self.y = [], []
            for i in range(n_spot):
                centroid, radius = self.get_centroid(self.found_spot[i])
                self.centroids.append(centroid)
                self.radiuses.append(radius)
                self.x.append(centroid.x)
                self.y.append(centroid.y)
        else :
            print(self.found_spot)
            centroid, radius = self.get_centroid(self.found_spot[0])
            self.centroids, self.radiuses = list(centroid), list(radius)
            self.x, self.y = list(centroid.x), list(centroid.y)
        self.get_substracted_background_image()
        positions = self.centroids
        apertures = [CircularAperture(position, r=radius) for position, radius in zip(positions, self.radiuses)]
        self.adu_count = [aperture.do_photometry(self.image)[0][0] for aperture in apertures] #No error propagation here
        return self.adu_count
        
    def get_detector_efficiency(self):
        data_dir = "/sdf/data/rubin/user/amouroux/DATA/CBP/config/"
        ETC = RubinCalibETC(self.exptime_config)
        ETC.get_detector_efficiency()
        d_eff = ETC.detector_efficiency
        w_index = np.where(ETC.rubin_wavelength == float(wavelength*1e9))[0][0]
        self.detector_efficiency = d_eff[w_index]

    def adu_to_electron(self, gain=1):
        self.el_count = self.adu_count*gain
        return self.el_count
        
    def electron_to_photon(self):
        self.h = 6.62607015e-34 #J s.
        self.c = 3e8 #m/s
        self.n_photons = self.el_count/self.detector_efficiency #photons
        return self.n_photons
    
    def get_flux(self, count):
        area = self.found_spot.getArea()*((1e-5)**2) #m^2
        self.flux = count/(area*total_exptime)
        return self.flux
        
    def get_physical_photon_flux(self):
        E_photon = self.h*self.c/self.wavelength # J
        E_tot = self.n_photons * self.E_photon
        #1px = .000010m
        area = self.found_spot.getArea()*((1e-5)**2) #m^2
        total_exptime = 16 #s #to correct for case dependency later
        self.flux = E_tot/(area*total_exptime) #W/m^2
        return get_physical_photon_flux