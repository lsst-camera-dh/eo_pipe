from collections import defaultdict
from typing import List

from astropy.table import Table as AstropyTable

from astro_metadata_translator import ObservationInfo

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT


class BiasStatsTaskConnections(pipeBase.PipelineTaskConnections,
                               dimensions=("instrument", "exposure",
                                           "detector")):
    raw_frames = cT.Input(name="raw_frame",
                          doc="Raw input frames",
                          storageClass="Exposure",
                          dimensions=("instrument", "exposure", "detector"),
                          multiple=True,
                          deferLoad=True)
    bias_stats = cT.Output(name="bias_stats",
                           doc="bias frame statistics",
                           storageClass="AstropyTable",
                           dimensions=("instrument", "exposure", "detector"))


class BiasStatsTaskConfig(pipeBase.PipelineTaskConfig,
                          pipelineConnections=BiasStatsTaskConnections):
    edge_buffer = pexConfig.Field(doc="Number of pixels to exclude around "
                                  "the edge of the serial overscan region.",
                                  dtype=int, default=2)


class BiasStatsTask(pipeBase.PipelineTask):

    ConfigClass = BiasStatsTaskConfig
    _DefaultName = "biasStatsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.edge_buffer = self.config.edge_buffer

    def run(self, exposures: List[afwImage.Exposure]) -> AstropyTable:
        data = defaultdict(list)
        for exposure in exposures:
            obs_info = ObservationInfo(exposure.getMetadata())
            det = exposure.getDetector()
            for amp in det:
                bbox = amp.getRawSerialOverscanBBox()
                bbox.grow(-self.edge_buffer)
                overscan = exposure.getImage()[bbox]
                flags = afwMath.MEDIAN | afwMath.STDEVCLIP
                stats = afwMath.makeStatistics(overscan, flags)
                data['run'].append(obs_info.science_program)
                data['exposure_id'].append(obs_info.exposure_id)
                data['det_name'].append(det.getName())
                data['amp'].append(amp.getName())
                data['median'].append(stats.getValue(afwMath.MEDIAN))
                data['stdevclip'].append(stats.getValue(afwMath.STDEVCLIP))
        return AstropyTable(data)
