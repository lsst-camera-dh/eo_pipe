import lsst.afw.math as afwMath
from lsst.afw.cameraGeom import utils as cgu
from lsst.obs.lsst import LsstCam
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT


class MyImageSource:
    isTrimmed = True
    background = 0.0
    def __init__(self, exposure_handles):
        self.exposure_handles = exposure_handles

    def getCcdImage(self, det, imageFactory, binSize=1, *args, **kwargs):
        ccdImage = self.exposure_handles[det.getId()].get().getImage()
        ccdImage = afwMath.binImage(ccdImage, binSize)
        return afwMath.rotateImageBy90(ccdImage,
                                       det.getOrientation().getNQuarter()), det


class DarkMosaicTaskConnections(pipeBase.PipelineTaskConnections,
                                dimensions=("instrument",)):
    dark_exposures = cT.Input(
        name="eoDarkISR",
        doc="Input pre-processed exposures to combine.",
        storageClass="Exposure",
        dimensions=("instrument", "exposure", "detector"),
        multiple=True,
        deferLoad=True)

    output_mosaic = cT.Output(
        name="eoDarkMosaic",
        doc="Full focal plane mosaic of dark frames with ISR applied.",
        storageClass="ImageF",
        dimensions=("instrument",))

class DarkMosaicTaskConfig(pipeBase.PipelineTaskConfig,
                           pipelineConnections=DarkMosaicTaskConnections):
    """Configuration for DarkMosaicTask."""
    binSize = pexConfig.Field(doc="Bin size in pixels for rebinning.",
                              dtype=int, default=4)

class DarkMosaicTask(pipeBase.PipelineTask):
    """Create mosaic of ISR'd dark frames for a particular exposure."""
    ConfigClass = DarkMosaicTaskConfig
    _defaultName = "darkMosaicTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.binSize = self.config.binSize

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        exposure_handles = {_.dataId['detector']: _ for _ in
                            butlerQC.get(inputRefs)['dark_exposures']}
        image = self.run(exposure_handles)
        butlerQC.put(image, outputRefs)

    def run(self, exposure_handles):
        camera = LsstCam.getCamera()
        image_source = MyImageSource(exposure_handles)
        detectorNameList = [camera[detector].getName() for detector in
                            exposure_handles]
        output_mosaic = cgu.showCamera(camera, imageSource=image_source,
                                       detectorNameList=detectorNameList,
                                       binSize=self.binSize)
        return pipeBase.Struct(output_mosaic=output_mosaic)
