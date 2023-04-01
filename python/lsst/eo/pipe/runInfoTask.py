from collections import defaultdict
import pandas as pd
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT


__all__ = ['RunInfoTask']


class RunInfoTaskConnections(pipeBase.PipelineTaskConnections,
                            dimensions=("instrument", "detector")):
    raws = cT.Input(name="raw",
                    doc="Raw input frames",
                    storageClass="Exposure",
                    dimensions=("instrument", "exposure", "detector"),
                    multiple=True,
                    deferLoad=True)

    run_info = cT.Output(name="run_info",
                         doc="eo_run_info",
                         storageClass="DataFrame",
                         dimensions=("instrument", "detector"))


class RunInfoTaskConfig(pipeBase.PipelineTaskConfig,
                        pipelineConnections=RunInfoTaskConnections):
    pass


class RunInfoTask(pipeBase.PipelineTask):
    ConfigClass = RunInfoTaskConfig
    _DefaultName = "runInfoTask"

    _exposure_fields = ['physical_filter', 'id', 'obs_id', 'exposure_time',
                        'dark_time', 'observation_type', 'observation_reason',
                        'day_obs', 'seq_num']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, raws):
        data = defaultdict(list)
        for handle in raws:
            data['detector'].append(handle.dataId['detector'])
            exp_dict = handle.dataId.records['exposure'].toDict()
            for field in self._exposure_fields:
                data[field].append(exp_dict[field])
            timespan = exp_dict['timespan']
            data['mjd_begin'].append(timespan.begin.mjd)
            data['mjd_end'].append(timespan.end.mjd)
        df = pd.DataFrame(data).sort_values(by=['detector', 'mjd_begin'])
        return pipeBase.Struct(run_info=df)
