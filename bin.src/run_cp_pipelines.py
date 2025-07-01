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
parser.add_argument('--bps_yaml', type=str, default=None,
                    help="Specific pipeline to run for a given run type")
parser.add_argument('--max_cpe_count', type=int, default=5,
                    help='Maximum number of CalledProcessErrors per pipeline')
parser.add_argument('--no_confirm', action='store_true', default=False,
                    help='Run pipelines without confirmation.')
args = parser.parse_args()

config_file = './cp_pipelines_config.yaml'
if not os.path.isfile(config_file):
    src = os.path.join(os.environ['EO_PIPE_DIR'], 'data', config_file)
    if os.environ['INSTRUMENT_NAME'] == "LATISS":
        src = src.replace("cp_pipelines_config", "cp_pipelines_latiss_config")
    shutil.copy(src, config_file)

cp_pipelines = CpPipelines(config_file, verbose=not args.laconic,
                           dry_run=args.dry_run, no_confirm=args.no_confirm,
                           max_cpe_count=args.max_cpe_count)
cp_pipelines.submit(args.run_type, bps_yaml=args.bps_yaml)
