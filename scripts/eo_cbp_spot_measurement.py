import lsst.daf.butler as daf_butler
from lsst.obs.lsst import LsstComCam #LsstCam
from lsst.eo.pipe.cbpspotMeasurementTask import SpotMeasurementTask, SpotMeasurementTaskConfig
from astropy.table import Table
import numpy as np
import sys, os
import importlib

### Here we hack through eo_pipe 

# Arguments, should be environment variables or command line arguments
os.environ["SPOT_RUN"] = "cbp_filter_scan"
os.environ["DAY_OBS"] = "20241210"
os.environ["OUTPUT"] = "u/amouroux/${PAYLOAD_NAME}"
os.environ["SPOT_RUN"] = "cbp_filter_scan"
os.environ["SPOT_RUN"] = "cbp_filter_scan"
os.environ["PHYSICAL_FILTER"] = "r"
os.environ["INSTRUMENT_NAME"] = "LSSTComCam"
os.environ["INSTRUMENT_CLASS_NAME"] = "LsstComCam"
os.environ["INSTRUMENT_CLASS"] = f"lsst.obs.lsst.{os.getenv("INSTRUMENT_CLASS_NAME")}"
os.environ["WEEKLY"] = "w_2025_13"
os.environ["EO_PIPE_INCOLLECTION"] = "u/tguillem/cbp_comcam/overscan_2D_w_2025_05_20250320b"
os.environ["PAYLOAD_MODIFIER"] = f"main_spots_{os.getenv("DAY_OBS")}"
os.environ["PAYLOAD_NAME"] = f"eo_cp_spot_measurements_{os.getenv("PAYLOAD_MODIFIER")}_{os.getenv("SPOT_RUN")}_{os.getenv("PHYSICAL_FILTER")}_{os.getenv("WEEKLY")}"
os.environ["OUTPUT"] = f"u/amouroux/{os.getenv("PAYLOAD_NAME")}"
os.environ["BUTLER_CONFIG"] = "/repo/main"
os.environ["DATAQUERY"] = f"instrument='{os.getenv("INSTRUMENT_NAME")}' and exposure.observation_reason='{os.getenv("SPOT_RUN")}' and band='{os.getenv("PHYSICAL_FILTER")}' and day_obs={os.getenv("DAY_OBS")}"

#intrument_module = importlib.import_module(os.getenv("INSTRUMENT_CLASS"))
#camera = intrument_module.getCamera()
camera = LsstComCam.getCamera()
aperture_list=[200] #px

## Set the config 
config = SpotMeasurementTaskConfig()
config.aperture = aperture_list[0]  # Set the default aperture to the first in the list
config.doForcedPhotometry = False

# Pass the config to the task
task = SpotMeasurementTask(config=config)

#aperture_list=[200] #px
#aperture_list = [50,100,150,200,300,400,500]
#exposure = [int(x.strip('[], ')) for x in exposure]
#print(exposure)
#print("forced photometry mode :", forced)
# Camera setup
#camera = LsstComCam.getCamera()

# Butler setup
#repo = "/repo/main"
#{collections = ['u/tguillem/cbp_comcam/overscan_2D_w_2025_05_20250320b'] #["LSSTComCam/nightlyValidation"]
butler = daf_butler.Butler(os.getenv("BUTLER_CONFIG"), collections=os.getenv("EO_PIPE_INCOLLECTION"))
#instrument='LSSTComCam'
where = os.getenv("DATAQUERY")
print(where)
refs = sorted(set(butler.registry.queryDatasets("postISRCCD", where=where)),
                  key=lambda ref: ref.dataId['exposure'])[:18]  # Limit to 100 exposures for testing
print("got datasetrefs", refs)
raws = [daf_butler.DeferredDatasetHandle(butler, ref, None) for ref in refs]

outputs = task.run(raws[:],camera)
base_path = f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/test_eo_task"
fits_path = f"{base_path}/{os.getenv("DAY_OBS")}_{aperture_list[0]}.fits"
print(outputs.spot_catalog)
# Write outputs to FITS
outputs.spot_catalog.writeFits(fits_path)

# Read and modify the table
#spot_table = Table.read(fits_path)
#print(f"Initial spot table: {spot_table}")

# Remove unnecessary columns
#spot_table.remove_columns(["coord_ra", "coord_dec", "parent", "footprint"])

# Overwrite the FITS file with the modified table
#spot_table.write(fits_path, overwrite=True)
#print(f"Modified spot table: {spot_table}")


#for exp in exposure :
    # Query datasets
#    where = ("instrument='LSSTComCam' and "
             #"detector=6 and "
#             f'exposure={exp}')
    
    
 #   print(f"Number of references found: {len(refs)}")
#    
    # Task setup
 #   task = SpotMeasurementTask()
    
  #  print(f"Raw datasets: {raws}")
    
    # File paths 
   # if refs[0].dataId['band'] == "white":
    #    band = "none"
    #else :
    #    band = refs[0].dataId['band']
    #base_path = f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/thibault_isr_{band}_spots"
    #if not os.path.exists(base_path):
    #    os.makedirs(base_path)
    #for aperture in aperture_list:
    #    if forced=="True":
    #        print("forced photometry mode")
    #        ref_spots_path = f"{base_path}/{day_obs}_ref_spots.fits"
     #       if not os.path.isfile(ref_spots_path):
      #          print("no ref spots file found.")
       #         break
        #    ref_spots = Table.read(ref_spots_path)
         #   ref_spots.sort("detector")
            #centroid = [ref_spots["x"], ref_spots["y"]]
          #  outputs = task.run(raws[:],camera,aperture, forced=True, ref_spots=ref_spots, repo=repo)
           # fits_path = f"{base_path}/{day_obs}_{exp}_{aperture}_forced_fitted_bkg.fits"
       # else : 
        #    outputs = task.run(raws[:],camera,aperture, forced=False, repo=repo)
        #    fits_path = f"{base_path}/{day_obs}_{exp}_{aperture}.fits"

        # Write outputs to FITS
        #outputs.spot_catalog.writeFits(fits_path)

        # Read and modify the table
        #spot_table = Table.read(fits_path)
        #print(f"Initial spot table: {spot_table}")

        # Remove unnecessary columns
        #spot_table.remove_columns(["coord_ra", "coord_dec", "parent", "footprint"])

        # Overwrite the FITS file with the modified table
        #spot_table.write(fits_path, overwrite=True)
        #print(f"Modified spot table: {spot_table}")

"""

for aperture in aperture_list:
        outputs = task.run(raws[:],camera,aperture)
        outputs.spot_catalog.writeFits(f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/optimise_aperture/{instrument}_spots_{day_obs}_{exposure}_{aperture}.fits")
        spot_table = Table.read(f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/optimise_aperture/{instrument}_spots_{day_obs}_{exposure}_{aperture}.fits")
        print(spot_table)
        spot_table.remove_columns(["coord_ra","coord_dec","parent", "footprint"])
        spot_table.write(f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/optimise_aperture/{instrument}_spots_{day_obs}_{exposure}_{aperture}.fits", overwrite = True)
        print(spot_table)"""