from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astro_metadata_translator import ObservationInfo

from lsst.afw.cameraGeom import utils as cgu
import lsst.afw.math as afwMath
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.daf.butler as daf_butler
import lsst.geom
from lsst.obs.lsst import LsstCam
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from .plotting import plot_focal_plane, append_acq_run
from .overscan_analysis import raft_oscan_correlations
from .dsref_utils import RaftOutputRefsMapper, get_plot_locations_by_dstype


__all__ = ['ReadNoiseTask', 'ReadNoiseFpPlotsTask', 'OverscanCorrelationsTask']


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


def get_plot_locations(repo, collections):
    dstypes = ('read_noise_plot', 'overscan_correlation_plot')
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
    dxy = pexConfig.Field(doc="Size in pixels of sub-regions.",
                          dtype=int, default=20)


class ReadNoiseTask(pipeBase.PipelineTask):

    ConfigClass = ReadNoiseTaskConfig
    _DefaultName = "readNoiseTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.edge_buffer = self.config.edge_buffer
        self.nsamp = self.config.nsamp
        self.dxy = self.config.dxy

    def run(self, raw_frames, ptc_results):
        data = defaultdict(list)
        for ds_handle in raw_frames:
            exposure = ds_handle.get()
            obs_info = ObservationInfo(exposure.getMetadata())
            det = exposure.getDetector()
            for amp in det:
                amp_name = amp.getName()
                gain = ptc_results.gain[amp_name]
                data['run'].append(obs_info.science_program)
                data['exposure_id'].append(obs_info.exposure_id)
                data['det_name'].append(det.getName())
                data['amp_name'].append(amp_name)

                bbox = amp.getRawSerialOverscanBBox()
                bbox.grow(-self.edge_buffer)
                overscan = exposure.getImage()[bbox]
                sampler = SubRegionSampler(overscan, self.dxy, self.dxy)
                stdevs = [
                    afwMath.makeStatistics(subregion, afwMath.STDEV).getValue()
                    for _, subregion in zip(range(self.nsamp), sampler)]
                data['read_noise'].append(gain*float(np.median(stdevs)))
                data['median'].append(gain*float(np.median(overscan.array)))
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
        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(111)
        title = append_acq_run(self, "Read Noise (e-)")
        plot_focal_plane(ax, amp_data, camera=camera, z_range=self.z_range,
                         scale_factor=self.z_scale_factor, title=title)
        return pipeBase.Struct(read_noise_plot=fig)


class OverscanCorrelationsTaskConnections(pipeBase.PipelineTaskConnections,
                                          dimensions=("instrument",)):
    raw_frames = cT.Input(
        name="raw",
        doc="Raw bias image.",
        storageClass="Exposure",
        dimensions=("instrument", "exposure", "detector"),
        multiple=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    overscan_correlation_plot = cT.Output(
        name="overscan_correlation_plot",
        doc="Raft-level correlations between amp overscan regions",
        storageClass="Plot",
        dimensions=("instrument", "detector"),
        multiple=True)


class OverscanCorrelationsTaskConfig(pipeBase.PipelineTaskConfig,
                                     pipelineConnections=OverscanCorrelationsTaskConnections):
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    oscan_buffer = pexConfig.Field(doc="Size of buffer region around the "
                                   "perimeter of the overscan to exclude "
                                   "from the cross-correlation calculations.",
                                   dtype=int, default=10)
    cmap = pexConfig.Field(doc="Matplotlib color map to use", dtype=str,
                           default='jet')


class OverscanCorrelationsTask(pipeBase.PipelineTask):
    ConfigClass = OverscanCorrelationsTaskConfig
    _DefaultName = "overscanCorrelationsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.xfigsize, self.config.yfigsize
        self.buffer = self.config.oscan_buffer
        self.cmap = self.config.cmap

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        camera = inputs['camera']
        raft_data = defaultdict(dict)
        raws = inputs['raw_frames']
        #
        # Get the exposure IDs, sort, and use the first one for the analysis.
        target = sorted(list(set(_.dataId['exposure'] for _ in raws)))[0]
        #
        # Organize refs by raft and slot.
        acq_run = None
        for ref in raws:
            if acq_run is None:
                acq_run = ref.dataId.records['exposure'].science_program
            if ref.dataId['exposure'] != target:
                continue
            det_name = ref.dataId.records['detector'].full_name
            if camera[det_name].getType() != cgu.DetectorType.SCIENCE:
                # Skip corner rafts.
                continue
            raft, slot = det_name.split('_')
            raft_data[raft][slot] = ref

        self.run(acq_run, target, raft_data, camera, butlerQC,
                 outputRefs.overscan_correlation_plot)

    def run(self, acq_run, exposure, raft_data, camera, butlerQC, output_refs):
        raft_output_refs_mapper = RaftOutputRefsMapper(camera)
        ref_map = raft_output_refs_mapper.create(output_refs)

        for raft, raws in raft_data.items():
            title = (f"Overscan correlations, {raft}, acq. run {acq_run}, "
                     f"{exposure}")
            fig, _ = raft_oscan_correlations(raws, camera, buffer=self.buffer,
                                             title=title, cmap=self.cmap,
                                             figsize=self.figsize)
            butlerQC.put(fig, ref_map[raft])
            plt.close()
