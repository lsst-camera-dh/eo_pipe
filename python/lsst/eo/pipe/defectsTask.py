from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd

import lsst.afw.detection as afw_detect
import lsst.afw.math as afw_math
from lsst.cp.pipe import MergeDefectsTaskConfig, MergeDefectsTask
import lsst.daf.butler as daf_butler
import lsst.geom
from lsst.ip.isr import Defects
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from .dh_utils import convert_amp_data_to_df
from .plotting import plot_focal_plane, hist_amp_data, append_acq_run
from .dsref_utils import get_plot_locations_by_dstype


__all__ = ["MergeSelectedDefectsTask", "BrightDefectsPlotsTask",
           "DarkDefectsPlotsTask", "VampireDefectsTask",
           "VampireDefectsPlotsTask"]


def get_amp_data(repo, collections):
    """Return amp-level pixel defects data."""
    butler = daf_butler.Butler(repo, collections=collections)

    amp_data = defaultdict(lambda: defaultdict(dict))
    defect_fields = {
        'bright_defects_results': ('bright_columns', 'bright_pixels'),
        'dark_defects_results': ('dark_columns', 'dark_pixels'),
        'vampire_defects_results': ('vampire_defects',)
    }
    for dstype, fields in defect_fields.items():
        dsrefs = list(set(butler.registry.queryDatasets(dstype,
                                                        findFirst=True)))
        if not dsrefs:
            continue
        df = butler.get(dsrefs[0])

        for _, row in df.iterrows():
            for field in fields:
                amp_data[field][row.det_name][row.amp_name] = row[field]
    return {field: dict(data) for field, data in amp_data.items()}


def get_plot_locations(repo, collections):
    dstypes = ('bright_columns_fp_plot', 'bright_pixels_fp_plot',
               'dark_columns_fp_plot', 'dark_pixels_fp_plot',
               'bright_columns_fp_hist', 'bright_pixels_fp_hist',
               'dark_columns_fp_hist', 'dark_pixels_fp_hist',
               'vampire_defects_fp_plot', 'vampire_defects_fp_hist')
    return get_plot_locations_by_dstype(repo, collections, dstypes)


def get_overlap_subcolumns(row, bbox):
    regions = []
    xmin = max(row.x0, bbox.minX)
    xmax = min(row.x0 + row.width - 1, bbox.maxX)
    ymin = max(row.y0, bbox.minY)
    ymax = min(row.y0 + row.height - 1, bbox.maxY)
    for x in range(xmin, xmax + 1):
        llc = lsst.geom.Point2I(x, ymin)
        urc = lsst.geom.Point2I(x, ymax)
        regions.append(lsst.geom.Box2I(llc, urc))
    return regions


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
            regions.extend(get_overlap_subcolumns(row, bbox))
        total_region_area += sum(_.area for _ in regions)

        bad_cols = set()
        for region in regions:
            if region.height > colthresh:
                bad_cols.add(region.minX)
        bad_columns[amp_name] = len(bad_cols)
        isolated_regions = [_ for _ in regions if _.minX not in bad_cols]
        # Exclude columns from bad pixel counts
        bad_pixels[amp_name] = sum(_.area for _ in isolated_regions)
    assert total_area == total_region_area
    return bad_columns, bad_pixels


class MergeSelectedDefectsTaskConfig(MergeDefectsTaskConfig):
    """Configuration for merging single exposure defects, with
    exposures selected by image_type.
    """
    imageType = pexConfig.Field(
        dtype=str,
        doc=("Image type of exposures to select the types of defects, "
             "i.e., image_type='dark' for bright defects, "
             "image_type='flat' for dark defects"),
        default=None)  # Make no selection on image type


class MergeSelectedDefectsTask(MergeDefectsTask):
    ConfigClass = MergeSelectedDefectsTaskConfig
    _DefaultName = "mergeSelectedDefects"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        camera = inputs["camera"]

        inputDefects = []
        for candidate in inputs["inputDefects"]:
            image_type = candidate.getMetadata()\
                                  .get('cpDefectGenImageType', 'UNKNOWN')
            if (self.config.imageType is None or
                image_type.lower() == self.config.imageType.lower()):
                inputDefects.append(candidate)

        if not inputDefects:
            outputs = pipeBase.Struct(mergedDefects=Defects())
        else:
            outputs = self.run(inputDefects, camera)
        butlerQC.put(outputs, outputRefs)


class DefectsPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                  dimensions=("instrument",)):
    defects = cT.Input(
        name="defects",
        doc="Defects found by cp_pipe",
        storageClass="Defects",
        dimensions=("instrument", "detector"),
        multiple=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",))


class DefectsPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                             pipelineConnections=DefectsPlotsTaskConnections):
    """Configuration for PixelDefectsPlotsTask"""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    acq_run = pexConfig.Field(doc="Acquisition run number.",
                              dtype=str, default="")


class DefectsPlotsTask(pipeBase.PipelineTask):
    """
    Base class for creating summary plots of pixel and column defects
    over the focal plane.
    """
    ConfigClass = DefectsPlotsTaskConfig
    _DefaultName = "defectsPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.yfigsize, self.config.xfigsize
        self.colthresh = self.config.colthresh

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        defects = {_.dataId['detector']: _ for _ in inputs['defects']}
        camera = inputs['camera']
        outputs = self.run(defects, camera)
        butlerQC.put(outputs, outputRefs)

    def run(self, defects, camera):
        return ()

    def make_defects_focalplane_plots(self, defects, camera, defect_type):
        outputs = {}
        column_data = {}
        pixel_data = {}
        for detector, handle in defects.items():
            defects = handle.get()
            det = camera[detector]
            det_name = det.getName()
            column_data[det_name], pixel_data[det_name] \
                = tabulate_defects(det, defects, colthresh=self.colthresh)
        outputs['column_data'] = column_data
        outputs['pixel_data'] = pixel_data
        xlabel = f"{defect_type} columns".lower()
        title = append_acq_run(self, xlabel.title())
        outputs['columns_fp_plot'], outputs['columns_fp_hist'] \
            = self._make_fp_plot(column_data, title, xlabel, camera)
        xlabel = f"{defect_type} pixels".lower()
        title = append_acq_run(self, xlabel.title())
        outputs['pixels_fp_plot'], outputs['pixels_fp_hist'] \
            = self._make_fp_plot(pixel_data, title, xlabel, camera)
        return outputs

    def _make_fp_plot(self, amp_data, title, xlabel, camera, z_range=None):
        mosaic = plt.figure(figsize=self.figsize)
        ax = mosaic.add_subplot(111)
        plot_focal_plane(ax, amp_data, camera=camera, title=title,
                         z_range=z_range, use_log10=True)
        hist = plt.figure()
        hist_amp_data(amp_data, xlabel, hist_range=z_range, use_log10=True,
                      title=title)
        return mosaic, hist


class BrightDefectsPlotsTaskConnections(DefectsPlotsTaskConnections):
    bright_columns_fp_plot = cT.Output(
        name="bright_columns_fp_plot",
        doc="Bright columns focal plane mosaic",
        storageClass="Plot",
        dimensions=("instrument",))

    bright_columns_fp_hist = cT.Output(
        name="bright_columns_fp_hist",
        doc="Bright columns focal plane histogram",
        storageClass="Plot",
        dimensions=("instrument",))

    bright_pixels_fp_plot = cT.Output(
        name="bright_pixels_fp_plot",
        doc="Bright pixels focal plane mosaic",
        storageClass="Plot",
        dimensions=("instrument",))

    bright_pixels_fp_hist = cT.Output(
        name="bright_pixels_fp_hist",
        doc="Bright pixels focal plane histogram",
        storageClass="Plot",
        dimensions=("instrument",))

    bright_defects_results = cT.Output(
        name="bright_defects_results",
        doc="Data frame of bright column and pixel defects results",
        storageClass="DataFrame",
        dimensions=("instrument",))


class BrightDefectsPlotsTaskConfig(
        DefectsPlotsTaskConfig,
        pipelineConnections=BrightDefectsPlotsTaskConnections):
    colthresh = pexConfig.Field(doc=("Minimum # continguous pixels "
                                     "defining a bright column."),
                                dtype=int, default=20)


class BrightDefectsPlotsTask(DefectsPlotsTask):
    """
    Class for creating summary plots of bright column and pixel defects.
    """
    ConfigClass = BrightDefectsPlotsTaskConfig
    _DefaultName = "brightDefectsPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, defects, camera):
        defect_type = "bright"
        outputs = self.make_defects_focalplane_plots(defects, camera,
                                                     defect_type)
        return pipeBase.Struct(
            bright_defects_results=convert_amp_data_to_df(
                {'bright_columns': outputs['column_data'],
                 'bright_pixels': outputs['pixel_data']}),
            bright_columns_fp_plot=outputs['columns_fp_plot'],
            bright_columns_fp_hist=outputs['columns_fp_hist'],
            bright_pixels_fp_plot=outputs['pixels_fp_plot'],
            bright_pixels_fp_hist=outputs['pixels_fp_hist'])


class DarkDefectsPlotsTaskConnections(DefectsPlotsTaskConnections):
    defects = cT.Input(
        name="defects",
        doc="Defects found by cp_pipe",
        storageClass="Defects",
        dimensions=("instrument", "detector", "physical_filter"),
        multiple=True,
        deferLoad=True)

    dark_columns_fp_plot = cT.Output(
        name="dark_columns_fp_plot",
        doc="Dark columns focal plane mosaic",
        storageClass="Plot",
        dimensions=("instrument",))

    dark_columns_fp_hist = cT.Output(
        name="dark_columns_fp_hist",
        doc="Dark columns focal plane histogram",
        storageClass="Plot",
        dimensions=("instrument",))

    dark_pixels_fp_plot = cT.Output(
        name="dark_pixels_fp_plot",
        doc="Dark pixels focal plane mosaic",
        storageClass="Plot",
        dimensions=("instrument",))

    dark_pixels_fp_hist = cT.Output(
        name="dark_pixels_fp_hist",
        doc="Dark pixels focal plane histogram",
        storageClass="Plot",
        dimensions=("instrument",))

    dark_defects_results = cT.Output(
        name="dark_defects_results",
        doc="Data frame of dark column and pixel defects results",
        storageClass="DataFrame",
        dimensions=("instrument",))


class DarkDefectsPlotsTaskConfig(
        DefectsPlotsTaskConfig,
        pipelineConnections=DarkDefectsPlotsTaskConnections):
    colthresh = pexConfig.Field(doc=("Minimum # continguous pixels "
                                     "defining a dark column."),
                                dtype=int, default=100)


class DarkDefectsPlotsTask(DefectsPlotsTask):
    """
    Class for creating summary plots of dark column and pixel defects.
    """
    ConfigClass = DarkDefectsPlotsTaskConfig
    _DefaultName = "darkDefectsPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, defects, camera):
        defect_type = "dark"
        outputs = self.make_defects_focalplane_plots(defects, camera,
                                                     defect_type)
        return pipeBase.Struct(
            dark_defects_results=convert_amp_data_to_df(
                {'dark_columns': outputs['column_data'],
                 'dark_pixels': outputs['pixel_data']}),
            dark_columns_fp_plot=outputs['columns_fp_plot'],
            dark_columns_fp_hist=outputs['columns_fp_hist'],
            dark_pixels_fp_plot=outputs['pixels_fp_plot'],
            dark_pixels_fp_hist=outputs['pixels_fp_hist'])


class VampireDefectsTaskConnections(
        pipeBase.PipelineTaskConnections,
        dimensions=("instrument", "detector", "physical_filter")):
    """Task to detect bright defects on combined flats."""
    flat = cT.Input(
        name="flat",
        doc="Combined flat constructed by cp_pipe",
        storageClass="Exposure",
        dimensions=("instrument", "detector", "physical_filter"),
        isCalibration=True
    )

    vampire_defects = cT.Output(
        name="vampire_defects",
        doc="Bright defects from combined flats",
        storageClass="Defects",
        dimensions=("instrument", "detector", "physical_filter"),
    )

    vampire_defect_catalog = cT.Output(
        name="vampire_defect_catalog",
        doc="Catalog of bright defects from combined flats",
        storageClass="DataFrame",
        dimensions=("instrument", "detector", "physical_filter"),
    )


class VampireDefectsTaskConfig(
        pipeBase.PipelineTaskConfig,
        pipelineConnections=VampireDefectsTaskConnections):
    background_size = pexConfig.Field(
        doc="Approximate size in pixels of cells used for background scaling.",
        dtype=int, default=64
    )
    threshold_type = pexConfig.ChoiceField(
        dtype=str,
        doc="Type of detection threshold.",
        default="VALUE",
        allowed={"VALUE": "Use fractional_threshold * clipped mean.",
                 "SIGMA": "Use nsigma * clipped stdev + clipped mean."}
    )
    fractional_threshold = pexConfig.Field(
        doc=("Fractional threshold above clipped mean for"
             "vampire pixel detection."),
        dtype=float, default=1.2
    )
    nsigma = pexConfig.Field(
        doc=("Number of clipped stdev above clipped mean for "
             "vampire pixel detection."),
        dtype=float, default=5.0
    )


class VampireDefectsTask(pipeBase.PipelineTask):
    """
    Class for detecting bright defects from combined flats.
    """
    ConfigClass = VampireDefectsTaskConfig
    _DefaultName = "vampireDefectsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.bg_size = self.config.background_size
        self.threshold_type = self.config.threshold_type
        self.frac_threshold = self.config.fractional_threshold
        self.nsigma = self.config.nsigma

    def run(self, flat):
        det = flat.getDetector()
        det_name = det.getName()
        image = flat.getMaskedImage()
        fp_list = []
        data = defaultdict(list)
        flags = afw_math.MEANCLIP | afw_math.STDEVCLIP
        for amp in det:
            amp_name = amp.getName()
            amp_image = self._background_scaling(image.Factory(image, amp.getBBox()))
            stats = afw_math.makeStatistics(amp_image, flags)
            mean = stats.getValue(afw_math.MEANCLIP)
            if self.threshold_type == "VALUE":
                threshold = afw_detect.createThreshold(self.frac_threshold*mean,
                                                       'value')
            elif self.threshold_type == "SIGMA":
                stdev = stats.getValue(afw_math.STDEVCLIP)
                threshold = afw_detect.createThreshold(mean + self.nsigma*stdev,
                                                       'value')
            footprints = afw_detect.FootprintSet(amp_image, threshold)\
                                   .getFootprints()
            fp_list.extend(footprints)
            # Subtract the baseline for flux calculations.
            amp_image -= mean
            for index, fp in enumerate(footprints):
                data['det_name'].append(det_name)
                data['amp_name'].append(amp_name)
                data['index'].append(index)
                centroid = fp.getCentroid()
                data['xpos'].append(centroid.x)
                data['ypos'].append(centroid.y)
                heavy_fp = afw_detect.HeavyFootprintF(fp, amp_image)
                data['flux'].append(sum(heavy_fp.getImageArray()))
                data['peak_flux'].append(max(heavy_fp.getImageArray()))
        return pipeBase.Struct(
            vampire_defects=Defects.fromFootprintList(fp_list),
            vampire_defect_catalog=pd.DataFrame(data)
        )

    def _background_scaling(self, amp_image):
        bbox = amp_image.getBBox()
        nx = bbox.getWidth() // self.bg_size
        ny = bbox.getHeight() // self.bg_size
        bg_ctrl = afw_math.BackgroundControl(nx, ny)
        bg_obj = afw_math.makeBackground(amp_image, bg_ctrl)
        amp_image /= bg_obj.getImageF("LINEAR")
        return amp_image


class VampireDefectsPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                         dimensions=("instrument",)):
    vampire_defect_catalogs = cT.Input(
        name="vampire_defect_catalogs",
        doc="Catalogs of bright defects from combined flats",
        storageClass="DataFrame",
        dimensions=("instrument", "detector", "physical_filter"),
        deferLoad=True,
        multiple=True
    )

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",))

    vampire_defects_fp_plot = cT.Output(
        name="vampire_defects_fp_plot",
        doc="Vampire defects focal plane mosaic",
        storageClass="Plot",
        dimensions=("instrument",))

    vampire_defects_fp_hist = cT.Output(
        name="vampire_defects_fp_hist",
        doc="Vampire defects focal plane histogram",
        storageClass="Plot",
        dimensions=("instrument",))


class VampireDefectsPlotsTaskConfig(
        pipeBase.PipelineTaskConfig,
        pipelineConnections=VampireDefectsPlotsTaskConnections):
    """Configuration for PixelDefectsPlotsTask"""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    zmin = pexConfig.Field(doc="Lower range of color bar for # defects per amp",
                           dtype=float, default=0)
    zmax = pexConfig.Field(doc="Upper range of color bar for # defects per amp",
                           dtype=float, default=100)
    acq_run = pexConfig.Field(doc="Acquisition run number.",
                              dtype=str, default="")


class VampireDefectsPlotsTask(pipeBase.PipelineTask):
    """
    Class for creating summary plots of vampire defects.
    """
    ConfigClass = VampireDefectsPlotsTaskConfig
    _DefaultName = "vampireDefectsPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.yfigsize, self.config.xfigsize
        self.z_range = self.config.zmin, self.config.zmax

    def run(self, vampire_defect_catalogs, camera):
        amp_data = defaultdict(dict)
        for handle in vampire_defect_catalogs:
            det_name = camera[handle.dataId['detector']].getName()
            df = handle.get()
            if len(df) == 0:
                continue
            for amp_name in set(df['amp_name']):
                amp_data[det_name][amp_name] \
                    = len(df.query(f"amp_name=='{amp_name}'"))

        mosaic = plt.figure(figsize=self.figsize)
        ax = mosaic.add_subplot(111)
        title = append_acq_run(self, "Bright defects from combined flats")
        plot_focal_plane(ax, amp_data, camera=camera, title=title,
                         z_range=self.z_range)
        hist = plt.figure()
        hist_amp_data(amp_data, "# defects", title=title)
        return pipeBase.Struct(
            vampire_defects_fp_plot=mosaic,
            vampire_defects_fp_hist=hist
        )
