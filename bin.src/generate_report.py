#!/usr/bin/env python

import os
import argparse
from lsst.eo.pipe import generate_report

parser = argparse.ArgumentParser()
parser.add_argument('--repo', help='data repository', type=str, default=None)
parser.add_argument('--acq_run', help='acquisition run', type=str)
parser.add_argument('--pattern', help='pattern to find collections', type=str,
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
acq_run = args.acq_run
pattern = args.pattern
staging_dir = args.staging_dir
htmldir = args.htmldir

if repo is None:
    repo = os.environ['BUTLER_CONFIG']

if args.payload_modifier is None:
    payload_modifier = os.environ.get("PAYLOAD_MODIFIER", "")
else:
    payload_modifier = args.payload_modifier

weekly = os.environ['WEEKLY']
if args.pattern is None:
    user = os.environ['USER']
    pattern = f"u/{user}/eo_*{payload_modifier}_{acq_run}_{weekly}"

generate_report(repo, pattern, acq_run, staging_dir=staging_dir,
                htmldir=htmldir, weekly=weekly+payload_modifier)
