#!/usr/bin/env python
import os
import argparse
import time
from lsst.eo.pipe import ingest_pd_data

parser = argparse.ArgumentParser()
parser.add_argument("run", type=str, help="acquisition run number")
parser.add_argument("--instrument", type=str, default="LSSTCam",
                    help="Instrument, e.g., LSSTCam, LSST-TS8")
parser.add_argument("--output_run", type=str, default=None,
                    help="Run collection to contain photodiode data")
parser.add_argument("--repo", type=str, default="embargo",
                    help="Data repository")
parser.add_argument("--target_repo", type=str, default="embargo",
                    help="Data repository")
parser.add_argument("--dry_run", action='store_true', default=False,
                    help="Print the ingest commands, but do not execute.")
parser.add_argument("--wait_time", type=float, default=5,
                    help="Minutes to wait between butler queries for pd data")

args = parser.parse_args()

while True:
    num_ingested = ingest_pd_data(args.run, instrument=args.instrument,
                                  output_run=args.output_run, repo=args.repo,
                                  target_repo=args.target_repo,
                                  dry_run=args.dry_run)
    if num_ingested == 0:
        break
    time.sleep(60*args.wait_time)

