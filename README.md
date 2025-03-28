This package contains the LSST Camera team code for electro-optical test code used for sensor testing.
The algorithms are based on the original LSST CCD testing specification, [LCA-128](https://docushare.lsst.org/docushare/dsweb/Get/LCA-128),
and the CCD Electro-Optical Testing Methods document [LCA-10103](https://docushare.lsst.org/docushare/dsweb/Get/LCA-10103).

The code in this package analyzes outputs from the LSST calibration pipelines, [cp_pipe](https://github.com/lsst/cp_pipe), and implements some custom pipelines using standard LSST code.
In particular, the instrument signature removal code is from the [ip_isr](https://github.com/lsst/ip_isr) package.  The cp_pipe and ip_isr repos are useful references for the default config settings.

This branch contains the tools to run the eo_task that will perform the photometry on the CBP images. 
1. scripts folder contains the python scripts to do the photometry
	1. eo_cbp_generate_ref_spots.py should be run first. It takes day_obs and band as keywords to find the exposures relative to a day of data acquisition in a band and finds the position of the spots on each detector. It then generates a table with the average position (x,y) of the 5 brightest spots found. The exposures are found via a table generated outside the module.
	2. eo_cp_spot_measurement_mproc.py should then be run. It takes day_obs, band and forced as keywords. Forced is likely to be set to "True", meaning that the photometry will be done in forced mode using the reference spots found by eo_cbp_generate_ref_spots.py. This code sends jobs in packets of five exposures, running the code from eo_cbp_spot_measurement.py
	3. eo_cbp_spot_measurement.py runs the eo task which performs the photometry using a list of exposures. It also transforms and saves the generated tables.
2. python/lsst/eo/pipe contains two Python scripts to perform the photometry and save the results.
	1. cbp_spotMeasurementTask.py performs the photometry by finding spots or using reference spots (if in forced mode) and generates a table containing the determined quantities such as spot position, total count, estimated mean background behind the spot, etc...
	2. photometry.py is the core package for performing the photometry. It is divided into three classes : Spot, ImageData and AperturePhotometry. Spot contains all the geometric information of the spot, while ImageData contains the processed image information (e.g  the image itself). AperturePhotometry likely takes an ImageData object and a Spot object as an input and perform the photometry with the "do_forced_aperture_photometry" method. Before integrating the photon count in the given aperture, this method fits a 2d background on the image and substracts the estimated background to the initial image.
