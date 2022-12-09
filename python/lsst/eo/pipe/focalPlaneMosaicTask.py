import lsst.afw.math as afwMath
from lsst.afw.cameraGeom import utils as cgu
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT


__all__ = ['FocalPlaneMosaicTask']


class ImageSource:
    isTrimmed = True
    background = 0.0

    def __init__(self, exposure_handles):
        self.exposure_handles = exposure_handles

    def getCcdImage(self, det, imageFactory, binSize=1, *args, **kwargs):
        ccdImage = self.exposure_handles[det.getId()].get().getImage()
        ccdImage = afwMath.binImage(ccdImage, binSize)
        return afwMath.rotateImageBy90(ccdImage,
                                       det.getOrientation().getNQuarter()), det


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
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

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


class FocalPlaneMosaicTask(pipeBase.PipelineTask):
    """Create mosaic of ISR'd CCD frames for a particular exposure."""
    ConfigClass = FocalPlaneMosaicTaskConfig
    _DefaultName = "focalPlaneMosaicTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.binSize = self.config.binSize

    def run(self, exposures, camera):
        exposure_handles = {}
        detectorNameList = []
        for handle in exposures:
            detector = handle.dataId['detector']
            exposure_handles[detector] = handle
            detectorNameList.append(camera[detector].getName())
        image_source = ImageSource(exposure_handles)
        output_mosaic = cgu.showCamera(camera, imageSource=image_source,
                                       detectorNameList=detectorNameList,
                                       binSize=self.binSize)
        return pipeBase.Struct(output_mosaic=output_mosaic)
