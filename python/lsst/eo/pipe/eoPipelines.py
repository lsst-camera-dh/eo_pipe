import os
import sys
import yaml
import subprocess
import click


__all__ = ['EoPipelines']


class EoPipelines:
    """
    Class to manage submission of eo_pipe pipelines of a particular
    run type, e.g., all the pipelines associated with a B-protocol or
    a PTC run.
    """
    def __init__(self, eo_pipeline_config, verbose=True, dry_run=False):
        """
        Parameters
        ----------
        eo_pipeline_config : str
            Configuration file containing the lists of pipelines and
            environment variables for the various run types.
        verbose : bool [True]
            Set to True for verbose output.
        dry_run : bool [False]
            Set to True to do a dry-run, listing pipelines to be run.
        """
        self.verbose = verbose
        self.dry_run = dry_run
        with open(eo_pipeline_config) as fobj:
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

    def submit(self, run_type):
        """
        Run bps submit for each pipeline of the specified run type.

        Parameters
        ----------
        run_type : str
            Type of run, e.g., b_protocol or ptc.
        """
        self._check_env_vars(run_type)
        eo_pipe_dir = os.environ['EO_PIPE_DIR']
        failed = []
        if self.verbose or self.dry_run:
            print("Running pipelines:")
            for pipeline in self.config[run_type]['pipelines']:
                print(f"  {pipeline}")
            print()
            self._print_inCollection()
        if self.dry_run:
            return
        if not click.confirm("Proceed?", default=True) or not self.verbose:
            print("Aborting runs.")
            return
        for pipeline in self.config[run_type]['pipelines']:
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
