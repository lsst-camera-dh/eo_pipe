from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astro_metadata_translator import ObservationInfo

import lsst.afw.math as afwMath
import lsst.daf.butler as daf_butler
import lsst.geom
from lsst.obs.lsst import LsstCam
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from .plotting import plot_focal_plane, hist_amp_data, append_acq_run
from .dsref_utils import get_plot_locations_by_dstype


__all__ = ['ReadNoiseTask', 'ReadNoiseFpPlotsTask']


def get_amp_data(repo, collections, camera=None):
    """Return amp-level read noise data"""
    if camera is None:
        camera = LsstCam.getCamera()
    butler = daf_butler.Butler(repo, collections=collections)
    dsrefs = list(set(butler.registry.queryDatasets('eo_read_noise',
                                                    findFirst=True)))
    amp_data = {'read_noise': defaultdict(dict),
                'bias_level': defaultdict(dict)}
    for dsref in dsrefs:
        det = camera[dsref.dataId['detector']]
        det_name = det.getName()
        for amp in det:
            amp_name = amp.getName()
            df = butler.get(dsref).query(f"amp_name=='{amp_name}'")
            for column in amp_data:
                try:
                    amp_data[column][det_name][amp_name] = np.median(df[column])
                except KeyError:
                    continue
    amp_data = {key: dict(value) for key, value in amp_data.items()}
    return amp_data


def get_plot_locations(repo, collections):
    dstypes = ('read_noise_plot', 'read_noise_hist')
    return get_plot_locations_by_dstype(repo, collections, dstypes)


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

    ptc_results = cT.Input(name="ptc_results",
                           doc="PTC fit results",
                           storageClass="PhotonTransferCurveDataset",
                           dimensions=("instrument", "detector"),
                           isCalibration=True)

    camera = cT.PrerequisiteInput(name="camera",
                                  doc="Camera used in observations",
                                  storageClass="Camera",
                                  isCalibration=True,
                                  dimensions=("instrument",))

    read_noise = cT.Output(name="eo_read_noise",
                           doc="read noise statistics",
                           storageClass="DataFrame",
                           dimensions=("instrument", "detector"))

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.usePtcGains:
            self.inputs.remove('ptc_results')


class ReadNoiseTaskConfig(pipeBase.PipelineTaskConfig,
                          pipelineConnections=ReadNoiseTaskConnections):
    edge_buffer = pexConfig.Field(doc="Number of pixels to exclude around "
                                  "the edge of the serial overscan region.",
                                  dtype=int, default=2)
    nsamp = pexConfig.Field(doc="Number of sub-regions to sample.",
                            dtype=int, default=100)
    dxy = pexConfig.Field(doc="Size in pixels of sub-regions.",
                          dtype=int, default=20)
    usePtcGains = pexConfig.Field(doc="Use gains from PTC dataset. If False, "
                                  "then use the gains from obs_lsst",
                                  dtype=bool, default=True)


class ReadNoiseTask(pipeBase.PipelineTask):

    ConfigClass = ReadNoiseTaskConfig
    _DefaultName = "readNoiseTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.edge_buffer = self.config.edge_buffer
        self.nsamp = self.config.nsamp
        self.dxy = self.config.dxy
        self.use_ptc_gains = self.config.usePtcGains

    def run(self, raw_frames, ptc_results=None, camera=None):
        data = defaultdict(list)
        for ds_handle in raw_frames:
            exposure = ds_handle.get()
            obs_info = ObservationInfo(exposure.getMetadata())
            det = exposure.getDetector()
            det_name = det.getName()
            for amp in det:
                amp_name = amp.getName()
                gain = camera[det_name][amp_name].getGain()
                if self.use_ptc_gains:
                    ptc_gain = ptc_results.gain[amp_name]
                    if not np.isnan(ptc_gain):
                        gain = ptc_gain
                data['run'].append(obs_info.science_program)
                data['exposure_id'].append(obs_info.exposure_id)
                data['det_name'].append(det_name)
                data['amp_name'].append(amp_name)

                bbox = amp.getRawSerialOverscanBBox()
                bbox.grow(-self.edge_buffer)
                overscan = exposure.getImage()[bbox]
                sampler = SubRegionSampler(overscan, self.dxy, self.dxy)
                stdevs = [
                    afwMath.makeStatistics(subregion, afwMath.STDEV).getValue()
                    for _, subregion in zip(range(self.nsamp), sampler)]
                data['read_noise'].append(gain*float(np.median(stdevs)))
                data['bias_level'].append(float(np.median(overscan.array)))
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
                                  dimensions=("instrument",))

    read_noise_plot = cT.Output(name="read_noise_plot",
                                doc="Focal plane plot of read noise values",
                                storageClass="Plot",
                                dimensions=("instrument",))

    read_noise_hist = cT.Output(name="read_noise_hist",
                                doc="Histogram of read noise values",
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
    acq_run = pexConfig.Field(doc="Acquisition run number.",
                              dtype=str, default="")


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
        read_noise_plot = plt.figure(figsize=self.figsize)
        ax = read_noise_plot.add_subplot(111)
        title = append_acq_run(self, "Read Noise (e-)")
        plot_focal_plane(ax, amp_data, camera=camera, z_range=self.z_range,
                         scale_factor=self.z_scale_factor, title=title)
        read_noise_hist = plt.figure()
        hist_amp_data(amp_data, "read noise (e-)", hist_range=self.z_range,
                      scale_factor=self.z_scale_factor, title=title)
        return pipeBase.Struct(read_noise_plot=read_noise_plot,
                               read_noise_hist=read_noise_hist)
