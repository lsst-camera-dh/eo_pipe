from collections import defaultdict
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astro_metadata_translator import ObservationInfo

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.daf.butler as daf_butler
import lsst.geom
from lsst.obs.lsst import LsstCam
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from lsst.eo.pipe.plotting import plot_focal_plane


__all__ = ['ReadNoiseTask', 'ReadNoiseFpPlotsTask']


def get_amp_data(repo, collections, camera=None):
    """Return amp-level read noise data"""
    if camera is None:
        camera = LsstCam.getCamera()
    butler = daf_butler.Butler(repo, collections=collections)
    dsrefs = list(set(butler.registry.queryDatasets('eo_read_noise',
                                                    findFirst=True)))
    amp_data = defaultdict(dict)
    for dsref in dsrefs:
        det = camera[dsref.dataId['detector']]
        det_name = det.getName()
        for amp in det:
            amp_name = amp.getName()
            df = butler.getDirect(dsref).query(f"amp_name=='{amp_name}'")
            amp_data[det_name][amp_name] = np.median(df['read_noise'])
    return {'read_noise': dict(amp_data)}


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

    read_noise = cT.Output(name="eo_read_noise",
                           doc="read noise statistics",
                           storageClass="DataFrame",
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
                data['amp_name'].append(amp.getName())

                bbox = amp.getRawSerialOverscanBBox()
                bbox.grow(-self.edge_buffer)
                overscan = exposure.getImage()[bbox]
                sampler = SubRegionSampler(overscan, self.dxy, self.dxy)
                stdevs = [
                    afwMath.makeStatistics(subregion, afwMath.STDEV).getValue()
                    for _, subregion in zip(range(self.nsamp), sampler)]
                data['read_noise'].append(float(np.median(stdevs)))
                data['median'].append(float(np.median(overscan.array)))
        return pipeBase.Struct(read_noise=pd.DataFrame(data))


class ReadNoiseFpPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                      dimensions=("instrument",)):
    read_noise_data = cT.Input(name="eo_read_noise",
                               doc="read noise statistics",
                               storageClass="DataFrame",
                               dimensions=("instrument", "detector"),
                               multiple=True,
                               deferLoad=True)

    camera = cT.PrerequisiteInput(name="camera",
                                  doc="Camera used in observations",
                                  storageClass="Camera",
                                  isCalibration=True,
                                  dimensions=("instrument",),
                                  lookupFunction=lookupStaticCalibration)

    read_noise_plot = cT.Output(name="read_noise_plot",
                                doc="Focal plane plot of read noise values",
                                storageClass="Plot",
                                dimensions=("instrument",))


class ReadNoiseFpPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                                 pipelineConnections=ReadNoiseFpPlotsTaskConnections):
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    zmin = pexConfig.Field(doc="Minimum of color bar range.",
                           dtype=float, default=0)
    zmax = pexConfig.Field(doc="Maximum of color bar range.",
                           dtype=float, default=20)
    zscale_factor = pexConfig.Field(doc=("Scale factor to apply to z-values."
                                         "This should be a float-convertable "
                                         "string so that formatting is taken "
                                         "care of automatically when "
                                         "rendering the label in matplotlib"),
                                    dtype=str, default="1")


class ReadNoiseFpPlotsTask(pipeBase.PipelineTask):
    """
    Create focal plane mosaics of the serial and parallel CTI results.
    """
    ConfigClass = ReadNoiseFpPlotsTaskConfig
    _DefaultName = 'readNoiseFpPlotsTask'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.yfigsize, self.config.xfigsize
        self.z_range = self.config.zmin, self.config.zmax
        self.z_scale_factor = self.config.zscale_factor

    def run(self, read_noise_data, camera):
        amp_data = defaultdict(dict)
        for ref in read_noise_data:
            det = camera[ref.dataId['detector']]
            det_name = det.getName()
            df0 = ref.get()
            for amp in det:
                amp_name = amp.getName()
                df = df0.query(f"amp_name == '{amp_name}'")
                amp_data[det_name][amp_name] = np.median(df['read_noise'])
        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(111)
        plot_focal_plane(ax, amp_data, camera=camera, z_range=self.z_range,
                         scale_factor=self.z_scale_factor,
                         title='Read Noise')
        return pipeBase.Struct(read_noise_plot=fig)