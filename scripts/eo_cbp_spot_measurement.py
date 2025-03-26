import lsst.daf.butler as daf_butler
from lsst.obs.lsst import LsstComCam #LsstCam
from lsst.eo.pipe.cbp_spotMeasurementTask import SpotMeasurementTask, SpotMeasurementTaskConfig
from astropy.table import Table
import numpy as np
import sys, os

# Arguments
exposure, day_obs, forced = sys.argv[1:-2], sys.argv[-2], sys.argv[-1]
aperture_list=[200] #px
#aperture_list = [50,100,150,200,300,400,500]
exposure = [int(x.strip('[], ')) for x in exposure]
print(exposure)
print("forced photometry mode :", forced)
# Camera setup
camera = LsstComCam.getCamera()

# Butler setup
repo = "/repo/main"
collections = ['u/tguillem/cbp_comcam/overscan_2D_w_2025_05_20250320b'] #["LSSTComCam/nightlyValidation"]
butler = daf_butler.Butler(repo, collections=collections)
instrument='LSSTComCam'

for exp in exposure :
    # Query datasets
    where = ("instrument='LSSTComCam' and "
             #"detector=6 and "
             f'exposure={exp}')
    
    refs = sorted(set(butler.registry.queryDatasets("postISRCCD", where=where)),
                  key=lambda ref: ref.dataId['exposure'])
    print(f"Number of references found: {len(refs)}")
    
    # Task setup
    task = SpotMeasurementTask()
    raws = [daf_butler.DeferredDatasetHandle(butler, ref, None) for ref in refs]
    print(f"Raw datasets: {raws}")
    
    # File paths 
    if refs[0].dataId['band'] == "white":
        band = "none"
    else :
        band = refs[0].dataId['band']
    base_path = f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/thibault_isr_{band}_spots"
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    for aperture in aperture_list:
        if forced=="True":
            print("forced photometry mode")
            ref_spots_path = f"{base_path}/{day_obs}_ref_spots.fits"
            if not os.path.isfile(ref_spots_path):
                print("no ref spots file found.")
                break
            ref_spots = Table.read(ref_spots_path)
            ref_spots.sort("detector")
            #centroid = [ref_spots["x"], ref_spots["y"]]
            outputs = task.run(raws[:],camera,aperture, forced=True, ref_spots=ref_spots, repo=repo)
            fits_path = f"{base_path}/{day_obs}_{exp}_{aperture}_forced_fitted_bkg.fits"
        else : 
            outputs = task.run(raws[:],camera,aperture, forced=False, repo=repo)
            fits_path = f"{base_path}/{day_obs}_{exp}_{aperture}.fits"

        # Write outputs to FITS
        outputs.spot_catalog.writeFits(fits_path)

        # Read and modify the table
        spot_table = Table.read(fits_path)
        print(f"Initial spot table: {spot_table}")

        # Remove unnecessary columns
        spot_table.remove_columns(["coord_ra", "coord_dec", "parent", "footprint"])

        # Overwrite the FITS file with the modified table
        spot_table.write(fits_path, overwrite=True)
        print(f"Modified spot table: {spot_table}")

"""

for aperture in aperture_list:
        outputs = task.run(raws[:],camera,aperture)
        outputs.spot_catalog.writeFits(f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/optimise_aperture/{instrument}_spots_{day_obs}_{exposure}_{aperture}.fits")
        spot_table = Table.read(f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/optimise_aperture/{instrument}_spots_{day_obs}_{exposure}_{aperture}.fits")
        print(spot_table)
        spot_table.remove_columns(["coord_ra","coord_dec","parent", "footprint"])
        spot_table.write(f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/optimise_aperture/{instrument}_spots_{day_obs}_{exposure}_{aperture}.fits", overwrite = True)
        print(spot_table)"""