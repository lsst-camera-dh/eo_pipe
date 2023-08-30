from lsst.afw.cameraGeom import utils as cgu
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from .plotting import ImageSource


__all__ = ['FocalPlaneMosaicTask']


class FocalPlaneMosaicTaskConnections(pipeBase.PipelineTaskConnections,
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
        dimensions=("instrument",))

    output_mosaic = cT.Output(
        name="eoFpMosaic",
        doc="Full focal plane mosaic of CCD exposures with ISR applied.",
        storageClass="ImageF",
        dimensions=("instrument",))


class FocalPlaneMosaicTaskConfig(pipeBase.PipelineTaskConfig,
                                 pipelineConnections=FocalPlaneMosaicTaskConnections):
    """Configuration for FocalPlaneMosaicTask."""
    binSize = pexConfig.Field(doc="Bin size in pixels for rebinning.",
                              dtype=int, default=4)
    filterWavelengthMap = pexConfig.DictField(
        doc=("Map of physical_filter name to wavelengths used for "
             "transmission curve integration"),
        keytype=str,
        itemtype=str,
        default={})
    use_ccob_led = pexConfig.Field(
        doc=("Use CCOBLED keyword instead of physical_filter for "
             "transmission curve scaling."),
        dtype=bool, default=True)


class FocalPlaneMosaicTask(pipeBase.PipelineTask):
    """Create mosaic of ISR'd CCD frames for a particular exposure."""
    ConfigClass = FocalPlaneMosaicTaskConfig
    _DefaultName = "focalPlaneMosaicTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.binSize = self.config.binSize
        self.filter_wl_map = self.config.filterWavelengthMap
        self.use_ccob_led = self.config.use_ccob_led

    def run(self, exposures, camera):
        exposure_handles = {}
        detectorNameList = []
        physical_filter = None
        for handle in exposures:
            if physical_filter is None:
                if self.use_ccob_led:
                    md = handle.butler.get(
                        handle.ref.makeComponentRef('metadata'))
                    physical_filter = md.get('CCOBLED', None)
                else:
                    physical_filter = handle.dataId['physical_filter']
            detector = handle.dataId['detector']
            exposure_handles[detector] = handle
            detectorNameList.append(camera[detector].getName())
        wls = eval(self.filter_wl_map.get(physical_filter, 'None'))
        image_source = ImageSource(exposure_handles, wls=wls)
        output_mosaic = cgu.showCamera(camera, imageSource=image_source,
                                       detectorNameList=detectorNameList,
                                       binSize=self.binSize)
        return pipeBase.Struct(output_mosaic=output_mosaic)
