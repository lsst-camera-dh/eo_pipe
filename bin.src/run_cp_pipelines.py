#!/usr/bin/env python
import os
import shutil
import argparse
from lsst.eo.pipe import CpPipelines

parser = argparse.ArgumentParser()
parser.add_argument('--run_type', default='b_protocol',
                    help='Type of analysis run',
                    choices=['b_protocol', 'ptc'])
parser.add_argument('--laconic', action='store_true', default=False,
                    help='Verbosity flag')
parser.add_argument('--dry_run', action='store_true', default=False,
                    help='Dry-run flag')
parser.add_argument('--max_cpe_count', type=int, default=5,
                    help='Maximum number of CalledProcessErrors per pipeline')
args = parser.parse_args()

config_file = './cp_pipelines_config.yaml'
if not os.path.isfile(config_file):
    src = os.path.join(os.environ['EO_PIPE_DIR'], 'data', config_file)
    if os.environ['INSTRUMENT_NAME'] == "LATISS":
        src = src.replace("cp_pipelines_config", "cp_pipelines_latiss_config")
    shutil.copy(src, config_file)

if args.run_type in ('b_protocol',):
    # Run pipelines sequentially
    sequential_pipelines = CpPipelines(config_file, verbose=not args.laconic,
                                       dry_run=args.dry_run,
                                       max_cpe_count=args.max_cpe_count)
    sequential_pipelines.submit(args.run_type)
elif args.run_type in ('ptc',):
    # Run pipelines concurrently
    concurrent_pipelines = EoPipelines(config_file, verbose=not args.laconic,
                                       dry_run=args.dry_run)
    concurrent_pipelines.submit(args.run_type)
