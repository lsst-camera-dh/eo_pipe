#!/usr/bin/env python
import os
import shutil
import argparse
from lsst.eo.pipe import EoPipelines

parser = argparse.ArgumentParser()
parser.add_argument('--run_type', default='b_protocol',
                    help='Type of analysis run',
                    choices=['b_protocol', 'ptc', 'ptc_subset',
                             'flat_gain_stability',
                             'dark_mosaic', 'run_info', 'read_noise_no_ptc'])
parser.add_argument('--laconic', action='store_true', default=False,
                    help='Verbosity flag')
parser.add_argument('--dry_run', action='store_true', default=False,
                    help='Dry-run flag')
parser.add_argument('--bps_yaml', type=str, default=None,
                    help='Specific pipeline to run for a given run type')
args = parser.parse_args()

config_file = './eo_pipe_config.yaml'
if not os.path.isfile(config_file):
    src = os.path.join(os.environ['EO_PIPE_DIR'], 'data', config_file)
    shutil.copy(src, config_file)

eo_pipelines = EoPipelines(config_file, verbose=not args.laconic,
                           dry_run=args.dry_run)

eo_pipelines.submit(args.run_type, bps_yaml=args.bps_yaml)
