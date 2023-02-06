from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from lsst.afw.cameraGeom import utils as cgu
from lsst.afw.math import flipImage
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from .focalPlaneMosaicTask import ImageSource
from lsst.eo.pipe.plotting import cmap_range


__all__ = ['RaftMosaicTask']


class RaftMosaicTaskConnections(pipeBase.PipelineTaskConnections,
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
    """Create raft mosaic of ISR'd CCD frames for a particular exposure."""
    ConfigClass = RaftMosaicTaskConfig
    _DefaultName = "raftMosaicTask"

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
        # Map the representative detector number to raft name.
        det_lists = defaultdict(list)
        for det in camera:
            det_name = det.getName()
            raft, _ = det_name.split('_')
            det_lists[raft].append(det.getId())
        detector_raft_map = {min(det_list): raft
                             for raft, det_list in det_lists.items()}

        # Loop over calib_type (i.e., "bias", "dark", "flat") and
        # the lists of calib exposure refs per detector.
        for calib_type, exp_refs in calibs.items():
            # Map output references for each raft.
            ref_map = {}
            for ref in output_refs[calib_type]:
                detector = ref.dataId['detector']
                if detector in detector_raft_map:
                    raft = detector_raft_map[detector]
                    ref_map[raft] = ref
            # Loop over rafts and create mosaics
            for raft, refs in exp_refs.items():
                image_source = ImageSource(refs)
                mosaic = cgu.showCamera(camera, imageSource=image_source,
                                        detectorNameList=list(refs.keys()),
                                        binSize=self.binSize)
                mosaic = flipImage(mosaic, flipLR=False, flipTB=True)
                imarr = mosaic.array
                raft_plot = plt.figure(figsize=self.figsize)
                ax = raft_plot.add_subplot(111)
                image = plt.imshow(imarr, interpolation='nearest',
                                   cmap=self.cmap)
                vmin, vmax = cmap_range(imarr, nsig=self.nsig)
                norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
                image.set_norm(norm)
                plt.title(f'{calib_type}, {raft}')
                plt.tick_params(axis='both', which='both', top=False,
                                bottom=False, left=False, right=False,
                                labelbottom=False, labelleft=False)

                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.1)
                plt.colorbar(image, cax=cax)

                butlerQC.put(raft_plot, ref_map[raft])
                plt.close()
