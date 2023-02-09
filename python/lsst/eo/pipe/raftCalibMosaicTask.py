from collections import defaultdict
import matplotlib.pyplot as plt
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from .dsref_utils import RaftOutputRefsMapper
from .plotting import make_mosaic


__all__ = ['RaftCalibMosaicTask']


class RaftCalibMosaicTaskConnections(pipeBase.PipelineTaskConnections,
                                     dimensions=("instrument",)):
    bias = cT.Input(
        name="bias",
        doc="Combined biases.",
        storageClass="Exposure",
        dimensions=("instrument", "detector"),
        isCalibration=True,
        multiple=True,
        deferLoad=True)

    dark = cT.Input(
        name="dark",
        doc="Combined darks.",
        storageClass="Exposure",
        dimensions=("instrument", "detector"),
        isCalibration=True,
        multiple=True,
        deferLoad=True)

    flat = cT.Input(
        name="flat",
        doc="Combined flats.",
        storageClass="Exposure",
        dimensions=("instrument", "physical_filter", "detector"),
        isCalibration=True,
        multiple=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    bias_mosaic = cT.Output(
        name="raft_bias_mosaic",
        doc="Raft-level bias mosaic",
        storageClass="Plot",
        dimensions=("instrument", "detector"),
        multiple=True)

    dark_mosaic = cT.Output(
        name="raft_dark_mosaic",
        doc="Raft-level dark mosaic",
        storageClass="Plot",
        dimensions=("instrument", "detector"),
        multiple=True)

    flat_mosaic = cT.Output(
        name="raft_flat_mosaic",
        doc="Raft-level flat mosaic",
        storageClass="Plot",
        dimensions=("instrument", "detector"),
        multiple=True)


class RaftCalibMosaicTaskConfig(pipeBase.PipelineTaskConfig,
                                pipelineConnections=RaftCalibMosaicTaskConnections):
    """Configuration for RaftCalibMosaicTask."""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    nsig = pexConfig.Field(doc="Number of stdevs below and above median "
                           "for color bar range", dtype=float, default=5.0)
    binSize = pexConfig.Field(doc="Bin size in pixels for rebinning.",
                              dtype=int, default=2)
    cmap = pexConfig.Field(doc="Matplotlib color map", dtype=str, default='hot')


class RaftCalibMosaicTask(pipeBase.PipelineTask):
    """Create raft mosaic of ISR'd CCD frames for a particular exposure."""
    ConfigClass = RaftCalibMosaicTaskConfig
    _DefaultName = "raftCalibMosaicTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.xfigsize, self.config.yfigsize
        self.nsig = self.config.nsig
        self.binSize = self.config.binSize
        self.cmap = self.config.cmap

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Gather input reference handles, organized by raft and detector,
        # and output refs over all detectors, keyed by calib_type.
        inputs = butlerQC.get(inputRefs)
        calibs = {}
        output_refs = {}
        for calib_type in ('bias', 'dark', 'flat'):
            # input reference handles
            calibs[calib_type] = defaultdict(dict)
            for ref in inputs[calib_type]:
                raft, detector = (ref.dataId.records['detector'].raft,
                                  ref.dataId['detector'])
                calibs[calib_type][raft][detector] = ref
            # output refs.  These are references per detector.
            output_refs[calib_type] = eval(f'outputRefs.{calib_type}_mosaic')

        camera = inputs['camera']

        # Do output ref persistence in the .run method since
        # collecting all of the raft-level plots in a single returned
        # pipeBase.Struct uses too much memory for the entire
        # focal plane.
        self.run(calibs, camera, butlerQC, output_refs)

    def run(self, calibs, camera, butlerQC, output_refs):
        raft_output_refs_mapper = RaftOutputRefsMapper(camera)
        # Loop over calib_type (i.e., "bias", "dark", "flat") and
        # the lists of calib exposure refs per detector.
        for calib_type, exp_refs in calibs.items():
            # Map output references for each raft.
            ref_map = raft_output_refs_mapper.create(output_refs[calib_type])
            # Loop over rafts and create mosaics
            for raft, refs in exp_refs.items():
                raft_plot = make_mosaic(refs, camera, self.binSize,
                                        self.figsize, self.cmap, self.nsig)
                plt.title(f'{calib_type}, {raft}')
                butlerQC.put(raft_plot, ref_map[raft])
                plt.close()
