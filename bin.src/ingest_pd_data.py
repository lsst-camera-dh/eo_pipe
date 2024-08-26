#!/usr/bin/env python
import os
import argparse
from lsst.eo.pipe import ingest_pd_data

parser = argparse.ArgumentParser()
parser.add_argument("run", type=str, help="acquisition run number")
parser.add_argument("--instrument", type=str, default="LSSTCam",
                    help="Instrument, e.g., LSSTCam, LSST-TS8")
parser.add_argument("--output_run", type=str, default=None,
                    help="Run collection to contain photodiode data")
parser.add_argument("--repo", type=str, default="/repo/embargo",
                    help="Data repository")
parser.add_argument("--s3_bucket", type=str, default="rubin-summit",
                    help="Data repository")

args = parser.parse_args()

ingest_pd_data(args.run, instrument=args.instrument,
               output_run=args.output_run, repo=args.repo,
               s3_bucket=args.s3_bucket)
