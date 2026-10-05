#!/usr/bin/env python

import os
import argparse
from lsst.eo.pipe import generate_report

parser = argparse.ArgumentParser()
parser.add_argument('--repo', help='data repository', type=str, default=None)
parser.add_argument('--dataset_label', help='dataset label', type=str,
                    default=None)
parser.add_argument('--pattern', help='pattern to find collections', type=str,
                    default=None)
parser.add_argument('--collections', help="list of collections", type=str,
                    default=None)
parser.add_argument('--staging_dir', help='staging directory for plots',
                    type=str, default='./eo_report_staging')
parser.add_argument('--htmldir', help='directory to contain html output',
                    type=str, default='/sdf/group/rubin/web_data/lsstcam')
parser.add_argument('--payload_modifier',
                    help=('string to append to the WEEKLY folder name '
                          'to indicate a modified payload'),
                    type=str, default=None)

args = parser.parse_args()

repo = args.repo
dataset_label = args.dataset_label
pattern = args.pattern
collections = args.collections.split(",") if args.collections is not None else None
staging_dir = args.staging_dir
htmldir = args.htmldir

if repo is None:
    repo = os.environ['BUTLER_CONFIG']

if args.payload_modifier is None:
    payload_modifier = os.environ.get("PAYLOAD_MODIFIER", "")
else:
    payload_modifier = args.payload_modifier

if dataset_label is None:
    dataset_label = os.environ['DATASET_LABEL']

weekly = os.environ['WEEKLY']
if args.pattern is None:
    user = os.environ['USER']
    pattern = f"u/{user}/eo_*{payload_modifier}_{dataset_label}_{weekly}"

generate_report(repo, pattern, dataset_label, staging_dir=staging_dir,
                htmldir=htmldir, weekly=weekly+payload_modifier,
                collections=collections)
