from collections import defaultdict

import matplotlib.pyplot as plt

from lsst.afw.cameraGeom import utils as cgu
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from .dsref_utils import RaftOutputRefsMapper, get_plot_locations_by_dstype
from .raft_level_correlations import raft_oscan_correlations, \
    raft_imaging_correlations


__all__ = ['ImagingCorrelationsTask', 'OverscanCorrelationsTask']


def get_plot_locations(repo, collections):
    dstypes = ('overscan_correlation_plot', 'imaging_correlation_plot')
    return get_plot_locations_by_dstype(repo, collections, dstypes)


class ImagingCorrelationsTaskConnections(pipeBase.PipelineTaskConnections,
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
        dimensions=("instrument",))

    imaging_correlation_plot = cT.Output(
        name="imaging_correlation_plot",
        doc="Raft-level correlations between amp imaging regions",
        storageClass="Plot",
        dimensions=("instrument", "detector"),
        multiple=True)


class ImagingCorrelationsTaskConfig(
        pipeBase.PipelineTaskConfig,
        pipelineConnections=ImagingCorrelationsTaskConnections):
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    edge_buffer = pexConfig.Field(doc="Size of buffer region around the "
                                  "perimeter of the imaging region to exclude "
                                  "from the cross-correlation calculations.",
                                  dtype=int, default=10)
    cmap = pexConfig.Field(doc="Matplotlib color map to use", dtype=str,
                           default='jet')


class ImagingCorrelationsTask(pipeBase.PipelineTask):
    ConfigClass = ImagingCorrelationsTaskConfig
    _DefaultName = "imagingCorrelationsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.xfigsize, self.config.yfigsize
        self.buffer = self.config.edge_buffer
        self.cmap = self.config.cmap

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        camera = inputs['camera']
        raft_data = defaultdict(lambda: defaultdict(dict))
        raws = [ref for ref in inputs['raw_frames']
                if ref.dataId.records['exposure'].observation_type=='flat']
        #
        # Get the exposure IDs, sort, and use the first two for the analysis.
        targets = sorted(list(set(_.dataId['exposure'] for _ in raws)))[:2]
        #
        # Organize refs by raft and slot.
        acq_run = None
        for ref in raws:
            if acq_run is None:
                acq_run = ref.dataId.records['exposure'].science_program
            exposure = ref.dataId['exposure']
            if exposure not in targets:
                continue
            det_name = ref.dataId.records['detector'].full_name
            if camera[det_name].getType() != cgu.DetectorType.SCIENCE:
                # Skip corner rafts.
                continue
            raft, slot = det_name.split('_')
            raft_data[raft][exposure][slot] = ref

        self.run(acq_run, targets, raft_data, camera, butlerQC,
                 outputRefs.imaging_correlation_plot)

    def run(self, acq_run, exposures, raft_data, camera, butlerQC, output_refs):
        raft_output_refs_mapper = RaftOutputRefsMapper(camera)
        ref_map = raft_output_refs_mapper.create(output_refs)

        for raft, raws in raft_data.items():
            title = (f"imaging region correlations, {raft}, "
                     f"acq. run {acq_run}, {exposures}")
            flat1_refs, flat2_refs = list(raws.values())
            fig, _ = raft_imaging_correlations(flat1_refs, flat2_refs, camera,
                                               buffer=self.buffer,
                                               title=title, cmap=self.cmap,
                                               figsize=self.figsize)
            butlerQC.put(fig, ref_map[raft])
            plt.close()


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
        dimensions=("instrument",))

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
        raws = [ref for ref in inputs['raw_frames']
                if ref.dataId.records['exposure'].observation_type=='bias']
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
