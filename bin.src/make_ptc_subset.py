#!/usr/bin/env python

import os
import sys
import argparse
import lsst.daf.butler as daf_butler
from lsst.eo.pipe import PtcSelection

parser = argparse.ArgumentParser()
parser.add_argument("run", type=str, help="PTC run to down-select")
parser.add_argument("--in_collection", type=str, default=None,
                    help="Input collection containing combined biases, "
                    "darks, and flats")
parser.add_argument("--repo", type=str, default="/repo/ir2",
                    help="Data repository")
parser.add_argument("--instrument", type=str, default="LSST-TS8",
                    help="Instrument, e.g., 'LSSTCam', 'LSST-TS8'")
parser.add_argument("--stride", type=int, default=5,
                    help="Stride to use in down-selection")
parser.add_argument("--flux_keyword", type=str, default="CCOBFLUX",
                    help="FITS header keyword to use to match flat pairs"
                    "and order the exposures for down-selection")
parser.add_argument("--detector", type=int, default=20,
                    help="detector to use to querying for exposures")

args = parser.parse_args()

if args.in_collection is None:
    user = os.environ["USER"]
    b_protocol_run = os.environ['B_PROTOCOL_RUN']
    weekly = os.environ['WEEKLY']
    in_collection = f"u/{user}/defects_{b_protocol_run}_{weekly}"
else:
    in_collection = args.in_collection

butler = daf_butler.Butler(args.repo)
try:
    butler.registry.queryCollections(in_collection)
except daf_butler.registry.MissingCollectionError:
    print(f"The default input collection {in_collection} does not exist.")
    print("You must provide a viable one using the --in_collection option.")
    sys.exit(0)

ptc_selection = PtcSelection(args.run, args.repo, args.instrument)
ptc_collection = ptc_selection.make_chained_ptc_collection(
    args.in_collection, stride=args.stride,
    flux_keyword=args.flux_keyword, detector=args.detector)

print("\nSet the following in your bash session:")
print(f"export PTC_CHAINED_COLLECTION={ptc_collection}")
