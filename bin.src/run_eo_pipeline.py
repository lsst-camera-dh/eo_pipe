#!/usr/bin/env python
import argparse
from lsst.eo.pipe import EoPipelines

parser = argparse.ArgumentParser()
parser.add_argument('--run_type', help='Type of analysis run',
                    default='b_protocol', type=str)
parser.add_argument('--eo_config', help='Run type config file',
                    default=None, type=str)
args = parser.parse_args()

config_file = args.eo_config
if config_file is None:
    config_file = os.path.join(os.environ['EO_PIPE_DIR'], 'data',
                               'eo_pipe_config.yaml')
eo_pipelines = EoPipelines(config_file)

eo_pipelines.submit(args.run_type)
