import os
import sys
import shutil
import time
import datetime
import yaml
import subprocess
import click


__all__ = ['EoPipelines', 'CpPipelines']


class PipelinesBase:
    """
    Class to manage submission of cp_pipe or eo_pipe pipelines of a
    particular run type, e.g., all the pipelines associated with a
    B-protocol or a PTC run.
    """
    def __init__(self, pipeline_config, verbose=True, dry_run=False,
                 log_dir='./logs'):
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
        log_idr : str ['./logs']
            Directory for log files that contain the `bps submit`
            screen output.
        """
        self.verbose = verbose
        self.dry_run = dry_run
        self.log_dir = log_dir
        os.makedirs(log_dir, exist_ok=True)
        with open(pipeline_config) as fobj:
            self.config = yaml.safe_load(fobj)
        self._copy_bps_folder()

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

    def _copy_bps_folder(self):
        local_bps_tree = './bps'
        if os.path.isdir(local_bps_tree):
            return
        main_bps_tree = os.path.join(os.environ['EO_PIPE_DIR'], 'bps')
        shutil.copytree(main_bps_tree, local_bps_tree)

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
            self._print_in_collection()
        if self.dry_run:
            return
        if not click.confirm("Proceed?", default=True) or not self.verbose:
            print("Aborting runs.")
            return
        self._run_pipelines(pipelines)

    def _print_in_collection(self):
        pass

    def _run_pipelines(self, pipelines):
        raise NotImplementedError

    def _bps_submit_command(self, bps_sub_folder, bps_yaml):
        root_dir = '.'
        log_file = os.path.join(self.log_dir, bps_yaml.replace('.yaml', '.log'))
        if os.path.isfile(log_file):
            os.remove(log_file)
        command = ' '.join(['bps', 'submit',
                            os.path.join(root_dir, bps_sub_folder, bps_yaml),
                            '2>&1 | tee', log_file])
        print('\n*****')
        print(command)
        print('*****')
        return command, log_file


class EoPipelines(PipelinesBase):
    """
    Class to manage bps submissions of eo_pipe pipelines.
    """
    def __init__(self, eo_pipeline_config, verbose=True, dry_run=False,
                 log_dir='./logs'):
        super().__init__(eo_pipeline_config, verbose=verbose, dry_run=dry_run,
                         log_dir=log_dir)

    def _print_in_collection(self):
        print("Using input collection:")
        command = (f"butler query-collections {os.environ['BUTLER_CONFIG']} "
                   f"{os.environ['EO_PIPE_INCOLLECTION']}")
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError:
            print("\nError querying for input collection. Aborting")
            sys.exit(1)
        print()

    def _run_pipelines(self, pipelines):
        failed = []
        for pipeline in pipelines:
            command, log_file = self._bps_submit_command('bps', pipeline)
            subprocess.check_call(command, shell=True)
            # Because the `bps submit` output is redirected through
            # `tee` in order to capture error messages in the log
            # file, bps failures do not return a non-zero return code.
            # Instead, we'll do `grep RuntimeError <log file>`, so
            # that the jobs for which that grep succeeds are the `bps
            # submit` jobs that failed.
            try:
                subprocess.check_call(["grep", "RuntimeError", log_file])
            except subprocess.CalledProcessError:
                pass
            else:
                failed.append(pipeline)
        if failed:
            print("Failed bps submissions:")
        for item in failed:
            print("  ", item)


class CpPipelines(PipelinesBase):
    """
    Class to manage sequential bps submission of cp_pipe pipelines.
    """
    def __init__(self, cp_pipeline_config, verbose=True, dry_run=False,
                 log_dir='./logs'):
        super().__init__(cp_pipeline_config, verbose=verbose, dry_run=dry_run)

    def _run_pipelines(self, pipelines):
        for pipeline in pipelines:
            command, _ = self._bps_submit_command('bps/cp_pipe', pipeline)
            output = subprocess.check_output(command, stderr=subprocess.STDOUT,
                                             shell=True, text=True).split("\n")
            for line in output:
                if line.startswith('Submit dir:'):
                    submit_dir = line.strip().split()[-1]
                if line.startswith('Run Id:'):
                    run_id = line.strip().split()[-1]
                    break
            print(f"Tracking {os.path.basename(pipeline)} run {run_id}")
            state = None
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
