from collections import defaultdict
import matplotlib.pyplot as plt
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from .dsref_utils import RaftOutputRefsMapper
from .plotting import make_mosaic


__all__ = ['RaftMosaicTask']


class RaftMosaicTaskConnections(pipeBase.PipelineTaskConnections,
                                dimensions=("instrument",)):
    exposures = cT.Input(
        name="postISRCCD",
        doc="Input pre-processed exposures to combine.",
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

    raft_mosaic_plot = cT.Output(
        name="eoRaftMosaic",
        doc="Raft-level mosaic of CCD exposures with ISR applied.",
        storageClass="ImageF",
        dimensions=("instrument", "detector"),
        multiple=True)


class RaftMosaicTaskConfig(pipeBase.PipelineTaskConfig,
                                 pipelineConnections=RaftMosaicTaskConnections):
    """Configuration for RaftMosaicTask."""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    nsig = pexConfig.Field(doc="Number of stdevs below and above median "
                           "for color bar range", dtype=float, default=5.0)
    binSize = pexConfig.Field(doc="Bin size in pixels for rebinning.",
                              dtype=int, default=2)
    cmap = pexConfig.Field(doc="Matplotlib color map", dtype=str, default='hot')


class RaftMosaicTask(pipeBase.PipelineTask):
    """Create mosaic of ISR'd CCD frames for a particular exposure."""
    ConfigClass = RaftMosaicTaskConfig
    _DefaultName = "raftMosaicTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.xfigsize, self.config.yfigsize
        self.nsig = self.config.nsig
        self.binSize = self.config.binSize
        self.cmap = self.config.cmap

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        camera = inputs['camera']
        raft_data = defaultdict(dict)
        for ref in inputs['exposures']:
            raft = ref.dataId.records['detector'].raft
            detector = ref.dataId['detector']
            raft_data[raft][detector] = ref

        self.run(raft_data, camera, butlerQC, outputRefs.raft_mosaic_plot)

    def run(self, raft_data, camera, butlerQC, output_refs):
        # Map the output references for each raft.
        raft_output_refs_mapper = RaftOutputRefsMapper(camera)
        ref_map = raft_output_refs_mapper.create(output_refs)

        for raft, exposure_refs in raft_data.items():
            raft_plot = make_mosaic(exposure_refs, camera, self.binSize,
                                    self.figsize, self.cmap, self.nsig)
            exposure = list(exposure_refs.values())[0].dataId['exposure']
            plt.title(f'{exposure}, {raft}')
            butlerQC.put(raft_plot, ref_map[raft])
            plt.close()
