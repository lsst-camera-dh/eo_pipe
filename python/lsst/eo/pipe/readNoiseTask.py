from collections import defaultdict
from typing import List

import numpy as np

from astro_metadata_translator import ObservationInfo

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT


__all__ = ['ReadNoiseTask', 'SubRegionSampler']


class SubRegionSampler:
    def __init__(self, image, dx, dy):
        self.image = image
        bbox = image.getBBox()
        self.x_range = bbox.getMinX(), bbox.getMaxX() - dx
        self.y_range = bbox.getMinY(), bbox.getMaxY() - dy
        self.extent = lsst.geom.Extent2I(dx, dy)

    def __next__(self):
        ix = int(np.random.randint(*self.x_range))
        iy = int(np.random.randint(*self.y_range))
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(ix, iy), self.extent)
        return self.image.Factory(self.image, bbox)

    def __iter__(self):
        return self


class ReadNoiseTaskConnections(pipeBase.PipelineTaskConnections,
                               dimensions=("instrument", "detector")):
    raw_frames = cT.Input(name="raw_frames",
                          doc="Raw input frames",
                          storageClass="Exposure",
                          dimensions=("instrument", "exposure", "detector"),
                          multiple=True,
                          deferLoad=True)
    read_noise = cT.Output(name="read_noise",
                           doc="read frame statistics",
                           storageClass="StructuredDataDict",
                           dimensions=("instrument", "detector"))


class ReadNoiseTaskConfig(pipeBase.PipelineTaskConfig,
                          pipelineConnections=ReadNoiseTaskConnections):
    edge_buffer = pexConfig.Field(doc="Number of pixels to exclude around "
                                  "the edge of the serial overscan region.",
                                  dtype=int, default=2)
    nsamp = pexConfig.Field(doc="Number of sub-regions to sample.",
                            dtype=int, default=100)
    dxy = pexConfig.Field(doc="Size in pixels of sub-regions",
                          dtype=int, default=20)


class ReadNoiseTask(pipeBase.PipelineTask):

    ConfigClass = ReadNoiseTaskConfig
    _DefaultName = "readNoiseTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.edge_buffer = self.config.edge_buffer
        self.nsamp = self.config.nsamp
        self.dxy = self.config.dxy

    def run(self, raw_frames: List[afwImage.Exposure]) -> pipeBase.Struct:
        data = defaultdict(list)
        for ds_handle in raw_frames:
            exposure = ds_handle.get()
            obs_info = ObservationInfo(exposure.getMetadata())
            det = exposure.getDetector()
            for amp in det:
                data['run'].append(obs_info.science_program)
                data['exposure_id'].append(obs_info.exposure_id)
                data['det_name'].append(det.getName())
                data['amp'].append(amp.getName())

                bbox = amp.getRawSerialOverscanBBox()
                bbox.grow(-self.edge_buffer)
                overscan = exposure.getImage()[bbox]
                sampler = SubRegionSampler(overscan, self.dxy, self.dxy)
                stdevs = [
                    afwMath.makeStatistics(subregion, afwMath.STDEV).getValue()
                    for _, subregion in zip(range(self.nsamp), sampler)]
                data['read_noise'].append(float(np.median(stdevs)))
                data['median'].append(float(np.median(overscan.array)))
        return pipeBase.Struct(read_noise=dict(data))
