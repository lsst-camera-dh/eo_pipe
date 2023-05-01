import os
import yaml
import subprocess


__all__ = ['EoPipelines']


class EoPipelines:
    """
    Class to manage submission of eo_pipe pipelines of a particular
    run type, e.g., all the pipelines associated with a B-protocol or
    a PTC run.
    """
    def __init__(self, eo_pipeline_config):
        """
        Parameters
        ----------
        eo_pipeline_config : str
            Configuration file containing the lists of pipelines and
            environment variables for the various run types.
        """
        with open(eo_pipeline_config) as fobj:
            self.config = yaml.safe_load(fobj)

    def _check_env_vars(self, run_type):
        # Check for required env vars for the specified run type.
        required = (self.config['baseline']['env_vars']
                    + self.config[run_type]['env_vars'])
        missing = [_ for _ in required if _ not in os.environ]
        if missing:
            raise RuntimeError("Missing required environment variables: "
                               f"{missing}")

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
        for pipeline in self.config[run_type]['pipelines']:
            command = ['bps', 'submit',
                       os.path.join(eo_pipe_dir, 'bps', pipeline)]
            print(' '.join(command))
            try:
                subprocess.check_call(command)
            except subprocess.CalledProcessError as eobj:
                failed.append((pipeline, eobj))
        for item in failed:
            print(item)
