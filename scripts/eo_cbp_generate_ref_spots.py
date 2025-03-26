import lsst.daf.butler as daf_butler
from lsst.obs.lsst import LsstComCam, LsstCam
from lsst.eo.pipe.cbp_spotMeasurementTask import SpotMeasurementTask, SpotMeasurementTaskConfig
from astropy.table import Table
import numpy as np
import sys, os
import time
import subprocess
from eo_cbp_spot_measurement_mproc import main

# Arguments
#exposure, day_obs = sys.argv[1:-1], sys.argv[-1]
#aperture_list=[200] #px
#aperture_list = [50,100,150,200,300,400,500]
#exposure = [int(x.strip('[], ')) for x in exposure]
#print(exposure)
# Camera setup
day_obs, band = sys.argv[1], sys.argv[2]
camera = LsstComCam.getCamera()

"""
def get_filter_bandpass(band):
    filter_path = "/sdf/group/rubin/user/amouroux/comissioning/cbp_analysis/Pointing & ghosting modeling notebook"
    data = np.loadtxt(f"{filter_path}/ideal_{band}.dat")
    wavelength = data[:, 0]  # Première colonne : Wavelength
    transmission = data[:, 1] 
    indexes = np.where(transmission>.99)
    filter_bandpass = [np.min(wavelength[indexes]), np.max(wavelength[indexes])]
    return filter_bandpass
"""

def produce_ref_spots(day_obs, band):
    print("produce_ref_spots")
    job_ids = main(day_obs, band, forced=False)
    wait_for_jobs_to_finish(job_ids)
    exp_tables_path = "/sdf/group/rubin/user/amouroux/comissioning/cbp_analysis/notebooks/comcam_analysis/exposures"
    spot_tables_path = "/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement"
    exp_table = Table.read(f"{exp_tables_path}/{band}_{day_obs}.fits")
    exposure_list = list(exp_table["exps"])
    spot_table={"signal":[],"detector":[],f"exposure":[],"x":[],"y":[]}
    for exp in exposure_list:
        spots = Table.read(f"{spot_tables_path}/{band}_spots/{day_obs}_{exp}_200.fits")
        for det in spots["det_name"]:
            spot_table["signal"].append(spots[spots["det_name"] == det]["signal"][0])
            spot_table["detector"].append(det)
            spot_table["exposure"].append(exp)
            spot_table["x"].append(spots[spots["det_name"] == det]["x"][0])
            spot_table["y"].append(spots[spots["det_name"] == det]["y"][0])
    spot_table = Table(spot_table)
    spot_table.sort("signal", reverse = True)
    brightest_spots = {"signal":[],"detector":[],f"exposure_list":[],"x":[],"y":[],"signal_list":[], "x_list":[],"y_list":[]}
    for det in spots["det_name"]:
        spot_table_det = spot_table[spot_table["detector"] == det][:5]
        brightest_spots["detector"].append(det)
        brightest_spots["exposure_list"].append(spot_table_det["exposure"].value)
        brightest_spots["signal_list"].append(spot_table_det["signal"].value)
        brightest_spots["x_list"].append(spot_table_det["x"].value)
        brightest_spots["y_list"].append(spot_table_det["y"].value)
        brightest_spots["signal"].append(np.mean(spot_table_det["signal"]))
        brightest_spots["x"].append(np.mean(spot_table_det["x"]))
        brightest_spots["y"].append(np.mean(spot_table_det["y"]))
    ref_path = f"{spot_tables_path}/raw_{band}_spots/{day_obs}_ref_spots.fits"
    brightest_spots = Table(brightest_spots)
    brightest_spots.write(ref_path, overwrite=True)

def wait_for_jobs_to_finish(job_ids, poll_interval=60):
    """
    Attend que les jobs SLURM spécifiés soient terminés.
    
    Args:
        job_ids (list): Liste des IDs des jobs SLURM à surveiller.
        poll_interval (int): Temps d'attente (en secondes) entre chaque vérification.
    """
    print(f"Attente de la fin des jobs : {', '.join(job_ids)}")
    
    while True:
        # Vérifier les jobs actifs avec `squeue`
        result = subprocess.run(["squeue", "--job", ",".join(job_ids)], stdout=subprocess.PIPE, text=True)
        squeue_output = result.stdout.strip()
        
        if len(squeue_output.splitlines()) <= 1:
            # S'il ne reste que la ligne d'en-tête, tous les jobs sont terminés
            print("Jobs done.")
            break
        else:
            # Afficher les jobs encore actifs (optionnel)
            print("Jobs running :")
            print(squeue_output)
        
        # Attendre avant la prochaine vérification
        time.sleep(poll_interval)

if __name__ == "__main__":
    produce_ref_spots(day_obs, band)