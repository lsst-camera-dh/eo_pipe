import subprocess
from astropy.table import Table
from pathlib import Path
import sys

def submit_job(exp, day, slurm_mem, slurm_output_path, current_dir, forced=False):
    """
    Submit a job to the Slurm workload manager.

    Parameters:
    exp (str): Exposure identifier.
    day (int): Observation day.
    slurm_mem (int): Memory allocation for the job in GB.
    slurm_output_path (str): Path to the Slurm output logs.
    current_dir (str): Current directory containing the script to run.
    """
    cmd = (
            f"sbatch --job-name=spot_measurement_mproc_{exp[0]} -t 0-30:00 -n 8 --mem {slurm_mem}G "
            f"-o {slurm_output_path}/{exp[0]}.out <<EOF\n"
            "#!/usr/bin/bash\n"
            f"source /sdf/group/rubin/user/amouroux/comissioning/setup.sh\n"
            f"setup eo_pipe -t amouroux\n"
            f"python {current_dir}/eo_cbp_spot_measurement.py {' '.join(map(str, exp))} {day} {forced}\n"
            "EOF"
        )
    try:
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        job_id = res.stdout.split("batch job")[1].split()[0]
        print(f"job_id:{job_id} submitted.  with forced = {forced}")
    except subprocess.CalledProcessError as e:
        print(f"Error submitting job for exposure {exp}: {e}")
    return job_id

def main(day_obs, filter_name, forced=False):
    slurm_mem = 20
    current_dir = "/sdf/home/a/amouroux/DATA/eo_pipe/scripts"
    exps_path = "/sdf/group/rubin/user/amouroux/comissioning/cbp_analysis/notebooks/comcam_analysis/exposures"
    exp_table = Table.read(f"{exps_path}/{filter_name}_{day_obs}.fits")
    exposures = list(exp_table[f"exposure"])
    print("Here are the exposures:", exposures)
    slurm_output_path = Path(f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/slurm_logs/{day_obs}/{filter_name}_spots/")
    slurm_output_path.mkdir(parents=True, exist_ok=True)
    job_id_list = []
    for i in range(int(len(exposures)/5)+1): # subdivide by packets of 5 exposures to not overload the slurm queue
        exps = exposures[i*5:(i+1)*5]
        #print(exps[0])
        if len(exps)==0:
            break
        print(exps,'\n',forced)
        job_id = submit_job(exps, day_obs, slurm_mem, slurm_output_path, current_dir, forced=forced)
        job_id_list.append(job_id)
    return job_id_list
    
if __name__ == "__main__":
    day_obs, filter_name, forced = sys.argv[1], sys.argv[2],  sys.argv[3]
    main(day_obs, filter_name, forced=forced)

"""slurm_mem = 10
current_dir = "/sdf/home/a/amouroux/DATA/eo_pipe/scripts"
not_processed = Table.read("/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/optimise_aperture/not_processed.fits")
for exp,aperture in zip(not_processed["exposure"],not_processed["aperture"]):
    slurm_output_path = f"/sdf/group/rubin/user/amouroux/DATA/CBP/Comcamspotmeasurement/slurm_logs/{str(exp)[:8]}/various_apertures/"
    cmd = f"sbatch --job-name=spot_measurement_mproc_{exp} -t  0-5:00 -n 2 --mem {slurm_mem}G -o {slurm_output_path}{exp}{aperture}.out <<?\n"
    cmd += "#!/usr/bin/bash\n"
    cmd += f"source /sdf/group/rubin/user/amouroux/comissioning/setup.sh\n"
    cmd += f"setup eo_pipe -t amouroux\n"
    cmd += f"python {current_dir}/eo_cbp_spot_measurement.py {exp} {str(exp)[:8]} {aperture}\n"
    res = subprocess.run(
    cmd, shell=True, capture_output=True, text=True, check=True
    )
    job_id = str(res.stdout).split("batch job")[1].split("\\")[0]
    print(f"job_id:{job_id} submitted.")"""