#!/usr/bin/env python
import argparse
import numpy as np
import lsst.daf.butler as daf_butler

parser = argparse.ArgumentParser()
parser.add_argument('run', type=str, help='acquisition run')
parser.add_argument('--exp_md_key', type=str, help='exposure metadata key',
                    default='science_program')
parser.add_argument('--repo', type=str, default='/repo/embargo',
                    help='data repository')
parser.add_argument('--instrument',
                    choices=['LSST-TS8', 'LSSTCam', 'LSSTComCam', 'LATISS'],
                    default='LSSTCam')

args = parser.parse_args()
instrument = args.instrument
butler = daf_butler.Butler(args.repo)

acq_run = args.run
exp_md_key = args.exp_md_key
if exp_md_key == "science_program":
    where = f"instrument='{instrument}' and exposure.{exp_md_key}='{acq_run}'"
else:
    where = f"instrument='{instrument}' and exposure.{exp_md_key}={acq_run}"
print(where)

refs = list(set(butler.registry.queryDatasets(
    "raw", where=where, collections=[f'{instrument}/raw/all']).expanded()))

refs = sorted(refs, key=lambda ref: ref.dataId['exposure'])
print(len(refs))
frame_data = []
exposure_rec = None
ptc_flats = set()
flats = set()
for ref in refs:
    exposure_rec = ref.dataId.records['exposure']
    md_key = str(eval(f"exposure_rec.{exp_md_key}"))
    frame_data.append((
        f"{md_key:6s}  "
        f"{ref.dataId['exposure']}  "
        f"{exposure_rec.obs_id:6s}  "
        f"{exposure_rec.observation_type:6s}  "
        f"{exposure_rec.observation_reason:6s}  "
        f"{ref.dataId['physical_filter']}  "))
    if (exposure_rec.observation_type == "flat" and
        exposure_rec.observation_reason == "flat"):
        ptc_flats.add(ref.dataId['exposure'])
    if exposure_rec.observation_type == "flat":
        flats.add(ref.dataId['exposure'])

unique_frames = np.unique(frame_data)
for i, item in enumerate(unique_frames):
    print(f"{i:05d}  {item}", flush=True)

if exposure_rec is not None:
    print(exposure_rec)

print()
print(len(set(butler.registry.queryDatasets("photodiode", where=where,
                                            collections=[f'{instrument}/photodiode']))),
      "photodiode datasets")
print(len(flats), "flat frames")
print(len(ptc_flats), "flat/flat frames")
