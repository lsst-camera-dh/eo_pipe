#!/usr/bin/env python
from collections import defaultdict
import argparse
import numpy as np
import lsst.daf.butler as daf_butler

parser = argparse.ArgumentParser()
parser.add_argument('run', type=str, help='acquisition run')
parser.add_argument('--exp_md_key', type=str, help='exposure metadata key',
                    default='science_program')
parser.add_argument('--repo', type=str, default='/repo/main',
                    help='data repository')
parser.add_argument('--instrument',
                    choices=['LSST-TS8', 'LSSTCam', 'LSSTComCam', 'LATISS'],
                    default='LSSTCam')
parser.add_argument('--get_ccobled', action="store_true", default=False,
                    help='include CCOBLED info for each exposure')

args = parser.parse_args()
instrument = args.instrument
get_ccobled = args.get_ccobled
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
high_flux_filter = None
low_flux_filter = None
ccd_count = defaultdict(lambda: 0)
for ref in refs:
    exposure = ref.dataId['exposure']
    ccd_count[exposure] += 1
    if ccd_count[exposure] > 1:
        continue
    if get_ccobled:
        md = butler.get(ref.makeComponentRef("metadata"))
        ccobled = md.get('CCOBLED', None)
    exposure_rec = ref.dataId.records['exposure']
    md_key = str(eval(f"exposure_rec.{exp_md_key}"))
    physical_filter = ref.dataId['physical_filter']
    entry = (
        f"{md_key:6s}  "
        f"{exposure}  "
        f"{exposure_rec.obs_id:6s}  "
        f"{exposure_rec.observation_type:6s}  "
        f"{exposure_rec.observation_reason:6s}  "
        f"{physical_filter}  ")
    if get_ccobled:
        entry += f"{ccobled}  "
    frame_data.append(entry)
    if (exposure_rec.observation_type == "flat" and
        exposure_rec.observation_reason == "flat"):
        ptc_flats.add(exposure)
    if exposure_rec.observation_type == "flat":
        flats.add(exposure)
    if exposure_rec.observation_reason == "sflat_hi":
        high_flux_filter = physical_filter
    if exposure_rec.observation_reason == "sflat_lo":
        low_flux_filter = physical_filter

unique_frames = np.unique(frame_data)
for i, item in enumerate(unique_frames):
    exposure = int(item.split()[1])
    print(f"{i:05d}  {item}", end="  ")
    print(f"{ccd_count[exposure]:3d}", flush=True)

if exposure_rec is not None:
    print()
    print("Example exposure record:")
    print(exposure_rec)

if high_flux_filter is not None and low_flux_filter is not None:
    print()
    print("For this run, these environment variables must be set:")
    print(f"export HIGH_FLUX_FILTER={high_flux_filter}")
    print(f"export LOW_FLUX_FILTER={low_flux_filter}")

print()
try:
    print(len(set(butler.registry.queryDatasets("photodiode", where=where,
                                                collections=[f'{instrument}/photodiode']))),
          "photodiode datasets")
except daf_butler.MissingCollectionError:
    print(f"No collection {instrument}/photodiode found.")

print(len(flats), "flat frames")
print(len(ptc_flats), "flat/flat frames")
print(len(unique_frames), "total frames")
