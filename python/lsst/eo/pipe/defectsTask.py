from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd

from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
from lsst.cp.pipe import MergeDefectsTaskConfig, MergeDefectsTask
import lsst.geom
from lsst.ip.isr import Defects
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from lsst.eo.pipe.plotting import plot_focal_plane


__all__ = ["MergeBrightDefectsTask", "MergeDarkDefectsTask", "DefectsPlotsTask"]


def get_overlap_region(row, bbox):
    llc = lsst.geom.Point2I(max(row.x0, bbox.minX),
                            max(row.y0, bbox.minY))
    urc = lsst.geom.Point2I(min(row.x0 + row.width - 1, bbox.maxX),
                            min(row.y0 + row.height - 1, bbox.maxY))
    return lsst.geom.Box2I(llc, urc)


def tabulate_defects(det, defects, colthresh=100):
    df0 = pd.DataFrame({_: defects.toDict()[_] for _ in
                        ('x0', 'y0', 'width', 'height')})
    total_area = sum(row.width*row.height for _, row in df0.iterrows())
    total_region_area = 0
    bad_columns = {}
    bad_pixels = {}
    for i, amp in enumerate(det):
        amp_name = amp.getName()
        bbox = amp.getBBox()
        df = df0.query(f'{bbox.minX} - width < x0 <= {bbox.maxX} and '
                       f'{bbox.minY} - height < y0 <= {bbox.maxY}')
        regions = []
        for _, row in df.iterrows():
            regions.append(get_overlap_region(row, bbox))
        total_region_area += sum(_.area for _ in regions)

        bad_cols = set()
        for region in regions:
            if region.height > colthresh:
                bad_cols.add(region.minX)
        bad_columns[amp_name] = len(bad_cols)
        isolated_regions = [_ for _ in regions if _.minX not in bad_cols]
        bad_pixels[amp_name] = (sum(_.area for _ in isolated_regions)
                                + bbox.height*len(bad_cols))
    assert total_area == total_region_area
    return bad_columns, bad_pixels


def convert_amp_data_to_df(amp_data, colname):
    data = defaultdict(list)
    for det_name, channel_data in amp_data.items():
        for amp_name, value in channel_data.items():
            data['det_name'].append(det_name)
            data['amp_name'].append(amp_name)
            data[colname].append(value)
    return pd.DataFrame(data)


class MergeBrightDefectsTask(MergeDefectsTask):
    ConfigClass = MergeDefectsTaskConfig
    _defaultName = "mergeBrightDefects"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        camera = inputs["camera"]

        inputDefects = []
        for candidate in inputs["inputDefects"]:
            image_type = candidate.getMetadata()\
                                  .get('cpDefectGenImageType', 'UNKNOWN')
            if image_type.lower() == "dark":
                inputDefects.append(candidate)

        if not inputDefects:
            outputs = pipeBase.Struct(mergedDefects=Defects())
        else:
            outputs = self.run(inputDefects, camera)
        butlerQC.put(outputs, outputRefs)


class MergeDarkDefectsTask(MergeDefectsTask):
    ConfigClass = MergeDefectsTaskConfig
    _defaultName = "mergeDarkDefects"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        camera = inputs["camera"]

        inputDefects = []
        for candidate in inputs["inputDefects"]:
            image_type = candidate.getMetadata()\
                                  .get('cpDefectGenImageType', 'UNKNOWN')
            if image_type.lower() != "dark":
                inputDefects.append(candidate)

        if not inputDefects:
            outputs = pipeBase.Struct(mergedDefects=Defects())
        else:
            outputs = self.run(inputDefects, camera)
        butlerQC.put(outputs, outputRefs)


class DefectsPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                  dimensions=("instrument",)):
    bright_defects = cT.Input(
        name="bright_defects",
        doc="Merged bright defects",
        storageClass="Defects",
        dimensions=("instrument", "detector"),
        multiple=True,
        isCalibration=True,
        deferLoad=True)

    dark_defects = cT.Input(
        name="dark_defects",
        doc="Merged dark defects",
        storageClass="Defects",
        dimensions=("instrument", "detector"),
        multiple=True,
        isCalibration=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    bright_columns_fp_plot = cT.Output(
        name="bright_columns_fp_plot",
        doc="Bright columns focal plane mosaic",
        storageClass="Plot",
        dimensions=("instrument",))

    bright_pixels_fp_plot = cT.Output(
        name="bright_pixels_fp_plot",
        doc="Bright pixels focal plane mosaic",
        storageClass="Plot",
        dimensions=("instrument",))

    dark_columns_fp_plot = cT.Output(
        name="dark_columns_fp_plot",
        doc="Dark columns focal plane mosaic",
        storageClass="Plot",
        dimensions=("instrument",))

    dark_pixels_fp_plot = cT.Output(
        name="dark_pixels_fp_plot",
        doc="Dark pixels focal plane mosaic",
        storageClass="Plot",
        dimensions=("instrument",))

    defects_results = cT.Output(
        name="defects_results",
        doc="Data frame of column and pixel defects results",
        storageClass="DataFrame",
        dimensions=("instrument",))


class DefectsPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                             pipelineConnections=DefectsPlotsTaskConnections):
    """Configuration for PixelDefectsPlotsTask"""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    bright_colthresh = pexConfig.Field(doc=("Minimum # continguous pixels "
                                            "defining a bright column."),
                                       dtype=int, default=20)
    dark_colthresh = pexConfig.Field(doc=("Minimum # continguous pixels "
                                          "defining a dark column."),
                                     dtype=int, default=100)

class DefectsPlotsTask(pipeBase.PipelineTask):
    """
    Summary plots of bright and dark pixel and column defects
    over the focal plane.
    """
    ConfigClass = DefectsPlotsTaskConfig
    _defaultName = "defectsPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.yfigsize, self.config.xfigsize
        self.bright_colthresh = self.config.bright_colthresh
        self.dark_colthresh = self.config.dark_colthresh

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        bright_defects = {_.dataId['detector']: _ for
                          _ in inputs['bright_defects']}

        dark_defects = {_.dataId['detector']: _ for
                        _ in inputs['dark_defects']}

        camera = inputs['camera']

        outputs = self.run(bright_defects, dark_defects, camera)
        butlerQC.put(outputs, outputRefs)

    def run(self, bright_defects, dark_defects, camera):
        outputs = {}

        bright_column_data = {}
        bright_pixel_data = {}
        for detector, handle in bright_defects.items():
            defects = handle.get()
            det = camera[detector]
            det_name = det.getName()
            columns, pixels = tabulate_defects(det, defects,
                                               colthresh=self.bright_colthresh)
            bright_column_data[det_name] = columns
            bright_pixel_data[det_name] = pixels

        outputs['bright_columns_fp_plot'] \
            = self._make_fp_plot(bright_column_data, "bright columns, run ?")
        outputs['bright_pixels_fp_plot'] \
            = self._make_fp_plot(bright_pixel_data, "bright pixels, run ?")

        dark_column_data = {}
        dark_pixel_data = {}
        for detector, handle in dark_defects.items():
            defects = handle.get()
            det = camera[detector]
            det_name = det.getName()
            columns, pixels = tabulate_defects(det, defects,
                                               colthresh=self.dark_colthresh)
            dark_column_data[det_name] = columns
            dark_pixel_data[det_name] = pixels

        outputs['dark_columns_fp_plot'] \
            = self._make_fp_plot(dark_column_data, "dark columns, run ?")
        outputs['dark_pixels_fp_plot'] \
            = self._make_fp_plot(dark_pixel_data, "dark pixels, run ?")

        outputs['defects_results'] = pd.concat(
            [convert_amp_data_to_df(bright_column_data, 'bright_columns'),
             convert_amp_data_to_df(bright_pixel_data, 'bright_pixels'),
             convert_amp_data_to_df(dark_column_data, 'dark_columns'),
             convert_amp_data_to_df(dark_pixel_data, 'dark_pixels')])

        return pipeBase.Struct(**outputs)

    def _make_fp_plot(self, amp_data, title, z_range=None):
        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(111)
        plot_focal_plane(ax, amp_data, title=title, z_range=z_range,
                         use_log10=True)
        return fig
