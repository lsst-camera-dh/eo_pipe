from collections import defaultdict
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.daf.butler as daf_butler
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from lsst.eo.pipe.plotting import plot_focal_plane


__all__ = ['BFAnalysisTask', 'BFAnalysisFpPlotsTask']


def get_amp_data(repo, collections):
    """Get Brighter-Fatter results for each amp in the focal plane,"""
    butler = daf_butler.Butler(repo, collections=collections)
    dsrefs = list(set(butler.registry.queryDatasets('bf_stats',
                                                    findFirst=True)))
    amp_data = defaultdict(lambda : defaultdict(dict))
    fields = ['bf_xcorr', 'bf_ycorr', 'bf_mean', 'bf_slope_x',
              'bf_slope_x_err', 'bf_slope_y', 'bf_slope_y_err']
    for dsref in dsrefs:
        df = butler.getDirect(dsref)
        for _, row in df.iterrows():
            for field in fields:
                amp_data[field][row.det_name][row.amp_name] = row[field]
    return {field: dict(data) for field, data in amp_data.items()}


class BFAnalysisTaskConnections(pipeBase.PipelineTaskConnections,
                                dimensions=("instrument", "detector")):
    ptc = cT.Input(name="ptc",
                   doc="PTC fit results",
                   storageClass="PhotonTransferCurveDataset",
                   dimensions=("instrument", "detector"),
                   isCalibration=True)

    camera = cT.PrerequisiteInput(name="camera",
                                  doc="Camera used in observations",
                                  storageClass="Camera",
                                  isCalibration=True,
                                  dimensions=("instrument",),
                                  lookupFunction=lookupStaticCalibration)

    bf_stats = cT.Output(name="bf_stats",
                         doc="Brighter-Fatter covariance statistics",
                         storageClass="DataFrame",
                         dimensions=("instrument", "detector"))

    bf_covariance_plots = cT.Output(name="bf_covariance_plots",
                                    doc=("Plots of covariance matrix elements "
                                         "vs mean signal each amp in a CCD"),
                                    storageClass="Plot",
                                    dimensions=("instrument", "detector"))


class BFAnalysisTaskConfig(pipeBase.PipelineTaskConfig,
                           pipelineConnections=BFAnalysisTaskConnections):
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=8)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=8)
    meanidx = pexConfig.Field(doc=("Index of covariance array to use for "
                                   "xcorr and ycorr"),
                              dtype=int, default=0)


class BFAnalysisTask(pipeBase.PipelineTask):
    """
    Task to characterize the brighter-fatter effect using covariance
    matrices computed from PTC data.
    """
    ConfigClass = BFAnalysisTaskConfig
    _DefaultName = "bfAnalysisTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.yfigsize, self.config.xfigsize
        self.meanidx = self.config.meanidx

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        # Get the detector number from the inputRefs.  There should
        # just be one.
        inputs['detector'] = inputRefs.ptc.dataId['detector']
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, ptc, camera, detector):
        det = camera[detector]
        det_name = det.getName()
        data = defaultdict(list)
        # Loop over amps and fit slopes to covariances vs mean signal.
        index = {}  # Indexes of covariance entries with outliers removed
        for amp in det:
            amp_name = amp.getName()
            means = np.array(ptc.rawMeans[amp_name])
            cov = np.array(ptc.covariances[amp_name])
            if len(cov) != len(means):
                continue
            index[amp_name] = ([i for i, use in
                                enumerate(ptc.expIdMask[amp_name]) if use],)
            means = means[index[amp_name]]
            cov = cov[index[amp_name]]
            if cov.size != 0:
                xcorr = cov[:, 1, 0]
                try:
                    slope_x, _, _, _, slope_x_err = stats.linregress(means,
                                                                     xcorr)
                except ValueError:
                    slope_x, slope_x_err = 0, 0
                ycorr = cov[:, 0, 1]
                try:
                    slope_y, _, _, _, slope_y_err = stats.linregress(means,
                                                                     ycorr)
                except ValueError:
                    slope_y, slope_y_err = 0, 0
            else:
                means = np.zeros(self.meanidx + 1)
                xcorr = np.zeros(self.meanidx + 1)
                ycorr = np.zeros(self.meanidx + 1)
                slope_x, slope_x_err, slope_y, slope_y_err = 0, 0, 0, 0
            data['det_name'].append(det_name)
            data['amp_name'].append(amp_name)
            data['bf_xcorr'].append(xcorr[self.meanidx])
            data['bf_ycorr'].append(ycorr[self.meanidx])
            data['bf_mean'].append(means[self.meanidx])
            data['bf_slope_x'].append(slope_x)
            data['bf_slope_x_err'].append(slope_x_err)
            data['bf_slope_y'].append(slope_y)
            data['bf_slope_y_err'].append(slope_y_err)
        bf_stats = pd.DataFrame(data)

        bf_covariance_plots = plt.figure(figsize=self.figsize)
        for i, (x, y) in enumerate(((1, 0), (2, 0), (0, 1), (0, 2), (1, 1)), 1):
            plt.subplot(3, 2, i)
            for amp_name, means in ptc.rawMeans.items():
                cov = np.array(ptc.covariances[amp_name])
                if len(cov) != len(means):
                    continue
                means = np.array(means)[index[amp_name]]
                cov = cov[index[amp_name]]
                if cov.size == 0:
                    continue
                plt.plot(means, cov[:, x, y]/means, label=amp_name)
            plt.title(f"cov{x}{y}")
            plt.legend(fontsize='x-small', ncol=4)
        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.suptitle(f"{det_name}")

        return pipeBase.Struct(bf_stats=bf_stats,
                               bf_covariance_plots=bf_covariance_plots)


class BFAnalysisFpPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                       dimensions=("instrument",)):
    bf_stats = cT.Input(
        name="bf_stats",
        doc="Brighter-Fatter covariance statistics",
        storageClass="DataFrame",
        dimensions=("instrument", "detector"),
        multiple=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    bf_xcorr_plot = cT.Output(
        name="bf_xcorr_plot",
        doc="Focal plane mosaic of cov10 for each amp",
        storageClass="Plot",
        dimensions=("instrument",))

    bf_ycorr_plot = cT.Output(
        name="bf_ycorr_plot",
        doc="Focal plane mosaic of cov01 for each amp",
        storageClass="Plot",
        dimensions=("instrument",))


class BFAnalysisFpPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                                  pipelineConnections=BFAnalysisFpPlotsTaskConnections):
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    zmin = pexConfig.Field(doc="Minimum of color bar range.",
                           dtype=float, default=-5)
    zmax = pexConfig.Field(doc="Maximum of color bar range.",
                           dtype=float, default=5)
    zscale_factor = pexConfig.Field(doc=("Scale factor to apply to z-values."
                                         "This should be a float-convertable "
                                         "string so that formatting is taken "
                                         "care of automatically when "
                                         "rendering the label in matplotlib"),
                                    dtype=str, default="1")


class BFAnalysisFpPlotsTask(pipeBase.PipelineTask):
    """
    Create focal plane mosaics of bf_xcorr and bf_ycorr.
    """
    ConfigClass = BFAnalysisFpPlotsTaskConfig
    _DefaultName = 'bfAnalysisFpPlotsTask'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.yfigsize, self.config.xfigsize
        self.z_range = self.config.zmin, self.config.zmax
        self.zscale_factor = self.config.zscale_factor

    def run(self, bf_stats, camera):
        # Unpack the bf_[xy]corr data for plotting.
        xcorr_data = defaultdict(dict)
        ycorr_data = defaultdict(dict)
        for handle in bf_stats:
            det = camera[handle.dataId['detector']]
            det_name = det.getName()
            df = handle.get()
            for _, row in df.iterrows():
                amp_name = row['amp_name']
                xcorr_data[det_name][amp_name] = row['bf_xcorr']
                ycorr_data[det_name][amp_name] = row['bf_ycorr']

        bf_xcorr_plot = plt.figure(figsize=self.figsize)
        ax = bf_xcorr_plot.add_subplot(111)
        plot_focal_plane(ax, xcorr_data, camera=camera, z_range=self.z_range,
                         scale_factor=self.zscale_factor, title='bf_xcorr')

        bf_ycorr_plot = plt.figure(figsize=self.figsize)
        ax = bf_ycorr_plot.add_subplot(111)
        plot_focal_plane(ax, ycorr_data, camera=camera, z_range=self.z_range,
                         scale_factor=self.zscale_factor, title='bf_ycorr')

        return pipeBase.Struct(bf_xcorr_plot=bf_xcorr_plot,
                               bf_ycorr_plot=bf_ycorr_plot)
