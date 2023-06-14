import os
import sys
import time
import datetime
import yaml
import subprocess
import click


__all__ = ['EoPipelines', 'CpPipelines']


class EoPipelines:
    """
    Class to manage submission of eo_pipe pipelines of a particular
    run type, e.g., all the pipelines associated with a B-protocol or
    a PTC run.
    """
    def __init__(self, pipeline_config, verbose=True, dry_run=False):
        """
        Parameters
        ----------
        pipeline_config : str
            Configuration file containing the lists of pipelines and
            environment variables for the various run types.
        verbose : bool [True]
            Set to True for verbose output.
        dry_run : bool [False]
            Set to True to do a dry-run, listing pipelines to be run.
        """
        self.verbose = verbose
        self.dry_run = dry_run
        with open(pipeline_config) as fobj:
            self.config = yaml.safe_load(fobj)

    def _check_env_vars(self, run_type):
        # Check for required env vars for the specified run type.
        required = (self.config['baseline']['env_vars']
                    + self.config[run_type]['env_vars'])
        missing = [_ for _ in required if _ not in os.environ]
        if self.verbose:
            print("Using environment variables:")
            for key in required:
                print(f"  {key}: {os.environ[key]}")
            print()
        if missing:
            raise RuntimeError("Missing required environment variables: "
                               f"{missing}")

    def _print_inCollection(self):
        print("Using input collection:")
        command = (f"butler query-collections {os.environ['BUTLER_CONFIG']} "
                   f"{os.environ['EO_PIPE_INCOLLECTION']}")
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError:
            print("\nError querying for input collection. Aborting")
            sys.exit(1)
        print()

    def submit(self, run_type, bps_yaml=None):
        """
        Run bps submit for each pipeline of the specified run type.

        Parameters
        ----------
        run_type : str
            Type of run, e.g., b_protocol or ptc.
        """
        self._check_env_vars(run_type)

        pipelines = self.config[run_type]['pipelines']
        if bps_yaml is not None:
            if bps_yaml not in pipelines:
                raise RuntimeError(f"{bps_yaml} not in {run_type} pipelines.")
            pipelines = [bps_yaml]

        if self.verbose or self.dry_run:
            print("Running pipelines:")
            for pipeline in pipelines:
                print(f"  {pipeline}")
            print()
            self._print_inCollection()
        if self.dry_run:
            return
        if not click.confirm("Proceed?", default=True) or not self.verbose:
            print("Aborting runs.")
            return
        self._run_pipelines(pipelines)

    def _run_pipelines(self, pipelines):
        failed = []
        eo_pipe_dir = os.environ['EO_PIPE_DIR']
        for pipeline in pipelines:
            command = ['bps', 'submit',
                       os.path.join(eo_pipe_dir, 'bps', pipeline)]
            print('\n*****')
            print(' '.join(command))
            print('*****')
            try:
                subprocess.check_call(command)
            except subprocess.CalledProcessError as eobj:
                failed.append((pipeline, eobj))
        for item in failed:
            print(item)


class CpPipelines(EoPipelines):
    """
    Class to manage sequential bps submission of cp_pipe pipelines.
    """
    def __init__(self, cp_pipeline_config, verbose=True, dry_run=False):
        super().__init__(cp_pipeline_config, verbose=verbose, dry_run=dry_run)

    def _print_inCollection(self):
        pass

    def _run_pipelines(self, pipelines):
        eo_pipe_dir = os.environ['EO_PIPE_DIR']
        for pipeline in pipelines:
            log_file = pipeline.replace('.yaml', '.log')
            if os.path.isfile(log_file):
                os.remove(log_file)
            command = ' '.join(['bps', 'submit',
                                os.path.join(eo_pipe_dir, 'bps', 'cp_pipe',
                                             pipeline), '2>&1 | tee', log_file])
            print('\n*****')
            print(command)
            print('*****')
            output = subprocess.check_output(command, stderr=subprocess.STDOUT,
                                             shell=True, text=True).split("\n")
            for line in output:
                if line.startswith('Submit dir:'):
                    submit_dir = line.strip().split()[-1]
                if line.startswith('Run Id:'):
                    run_id = line.strip().split()[-1]
                    break
            print(f"Tracking {os.path.basename(pipeline)} run {run_id}")
            while True:
                i = 0
                command = ["bps", "report", "--id", submit_dir]
                output = subprocess.check_output(command,
                                                 stderr=subprocess.STDOUT,
                                                 text=True).split('\n')
                for line in output:
                    if run_id in line and not line.startswith('Global job id'):
                        state = line[3:].strip().split()[0]
                        now = datetime.datetime\
                                      .now().strftime('%Y-%m-%d %H:%M:%S')
                        print(now, state)
                time.sleep(10)
                i += 1
                if state == 'SUCCEEDED':
                    break
                if state == 'FAILED':
                    raise RuntimeError(f"Pipeline {os.path.basename(pipeline)} "
                                       "failed")
