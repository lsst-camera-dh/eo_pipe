from collections import defaultdict
from astropy.io import fits
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
    metadata_fields = pexConfig.ListField(
        doc="Keyword names of values to extract from raw exposure metadata.",
        dtype=str, default=[])
    use_astropy = pexConfig.Field(
        doc="Flag to use astropy.io.fits to read the FITS header "
        "instead of using the .getMetadata function.",
        dtype=bool, default=True)


class RunInfoTask(pipeBase.PipelineTask):
    ConfigClass = RunInfoTaskConfig
    _DefaultName = "runInfoTask"

    _exposure_fields = ['physical_filter', 'id', 'obs_id', 'exposure_time',
                        'dark_time', 'observation_type', 'observation_reason',
                        'day_obs', 'seq_num']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.keys = self.config.metadata_fields
        self.use_astropy = self.config.use_astropy

    def run(self, raws):
        data = defaultdict(list)
        self.log.info("Found %d raw files", len(raws))
        for i, handle in enumerate(raws):
            self.log.info("processing %d: %s", i, handle.dataId['exposure'])
            data['detector'].append(handle.dataId['detector'])
            exp_dict = handle.dataId.records['exposure'].toDict()
            for field in self._exposure_fields:
                data[field].append(exp_dict[field])
            timespan = exp_dict['timespan']
            data['mjd_begin'].append(timespan.begin.mjd)
            data['mjd_end'].append(timespan.end.mjd)
            if self.keys:
                if self.use_astropy:
                    filepath = handle.butler.getURI(handle.ref).path
                    with fits.open(filepath) as hdus:
                        md = dict(hdus[0].header)
                else:
                    md = handle.get().getMetadata()
                for key in self.keys:
                    data[key].append(md.get(key, None))
        df = pd.DataFrame(data).sort_values(by=['detector', 'mjd_begin'])
        return pipeBase.Struct(run_info=df)
