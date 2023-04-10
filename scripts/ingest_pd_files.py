"""
Example script for ingesting photodiode readings files.
"""

import os
import glob
import subprocess

day_obs_folders = ['/fs/ddn/sdf/group/lsst/camera/IandT/R_and_D/ts8/20230310',
                   '/fs/ddn/sdf/group/lsst/camera/IandT/R_and_D/ts8/20230321']

pd_files = []
for folder in day_obs_folders:
    pd_files.extend(glob.glob(os.path.join(folder, 'TS_C*', 'Photodiode*')))

output_repo = '/repo/ir2'
output_run = 'LSST-TS8/photodiode'

print(len(pd_files))
for i, pd_file in enumerate(pd_files):
    if i % 10:
        print('.', end='', flush=True)
    else:
        print(i, end='', flush=True)
    command = f"butler ingest-photodiode --output-run {output_run} {output_repo} LSST-TS8 {pd_file}"
    #print(command)
    subprocess.check_call(command, shell=True)
print()
