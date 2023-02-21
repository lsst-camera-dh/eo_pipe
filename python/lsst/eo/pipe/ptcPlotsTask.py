from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import lsst.afw.math as afwMath
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
from lsst.cp.pipe.utils import funcAstier, funcPolynomial
import lsst.daf.butler as daf_butler
from lsst.obs.lsst import LsstCam
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from lsst.eo.pipe.plotting import plot_focal_plane


__all__ = ['PtcPlotsTask', 'PtcFpPlotsTask', 'RowMeansVarianceTask',
           'RowMeansVarianceFpPlotTask']


def get_amp_data(repo, collections, camera=None):
    """Get PTC results for each amp in the focal plane."""
    if camera is None:
        camera = LsstCam.getCamera()

    butler = daf_butler.Butler(repo, collections=collections)

    amp_data = defaultdict(lambda: defaultdict(dict))

    # Extra PTC results.
    dsrefs = list(set(butler.registry.queryDatasets('ptc', findFirst=True)))
    for dsref in dsrefs:
        det = camera[dsref.dataId['detector']]
        det_name = det.getName()
        ptc = butler.getDirect(dsref)
        for amp_name, values in ptc.ptcFitPars.items():
            if len(values) != 3:
                continue
            ptc_a00, ptc_gain, ptc_var = values
            amp_data['ptc_a00'][det_name][amp_name] = -ptc_a00
            amp_data['ptc_gain'][det_name][amp_name] = ptc_gain
            if ptc_var > 0:
                amp_data['ptc_noise'][det_name][amp_name] = np.sqrt(ptc_var)
            amp_data['ptc_turnoff'][det_name][amp_name] \
                = (np.max(ptc.finalMeans[amp_name])
                   if ptc.finalMeans[amp_name] else -1)

    # Extract row means variance slopes
    dsrefs = list(set(butler.registry.queryDatasets('row_means_variance_stats')))
    for ref in dsrefs:
        df = butler.getDirect(ref)
        for _, row in df.iterrows():
            amp_data['row_mean_var_slope'][row.det_name][row.amp_name] \
                = row.slope

    return {field: dict(data) for field, data in amp_data.items()}


class PtcPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                              dimensions=("instrument", "detector")):
    ptc_results = cT.Input(name="ptc_results",
                           doc="PTC fit results",
                           storageClass="PhotonTransferCurveDataset",
                           dimensions=("instrument", "detector"),
                           isCalibration=True,
                           deferLoad=True)
    camera = cT.PrerequisiteInput(name="camera",
                                  doc="Camera used in observations",
                                  storageClass="Camera",
                                  isCalibration=True,
                                  dimensions=("instrument",),
                                  lookupFunction=lookupStaticCalibration)
    ptc_plots = cT.Output(name="ptc_plots",
                          doc=("Plots of photon transfer curves for each "
                               "amp in a CCD"),
                          storageClass="Plot",
                          dimensions=("instrument", "detector"))


class PtcPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                         pipelineConnections=PtcPlotsTaskConnections):
    """Configuration for PtcPlotsTask."""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)


_ptc_func = {'EXPAPPROXIMATION': funcAstier,
             'POLYNOMIAL': funcPolynomial}


class PtcPlotsTask(pipeBase.PipelineTask):
    """Plot the photon transfer curves for all of the amps in a CCD."""
    ConfigClass = PtcPlotsTaskConfig
    _DefaultName = "ptcPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = (self.config.yfigsize, self.config.xfigsize)

    def run(self, ptc_results, camera):
        detector = ptc_results.dataId['detector']
        det = camera[detector]
        ptc = ptc_results.get()

        fig = plt.figure(figsize=self.figsize)
        for i, amp in enumerate(det, 1):
            plt.subplot(4, 4, i)
            amp_name = amp.getName()
            index = ([i for i, use in enumerate(ptc.expIdMask[amp_name])
                      if use],)
            means = np.array(ptc.rawMeans[amp_name])
            variances = np.array(ptc.rawVars[amp_name])
            plt.scatter(means, variances, s=8, color='blue')
            plt.scatter(means[index], variances[index], s=2, color='red')
            if ptc.ptcFitType == "EXPAPPROXIMATION":
                a00, gain, _ = ptc.ptcFitPars[amp_name]
                a00_err, gain_err, _ = ptc.ptcFitParsError[amp_name]
                turnoff = ptc.ptcTurnoff[amp_name]
                note = (f"{amp_name}\n"
                        f"gain = {gain:.2f} +/- {gain_err:.2f}\n"
                        f"a00 = {a00:.1e} +/- {a00_err:.1e}\n"
                        f"turnoff = {turnoff:7.0f}")
                plt.annotate(note, (0.05, 0.95), xycoords='axes fraction',
                             verticalalignment='top', size='x-small')
            if means.size != 0:
                if ptc.ptcFitType in _ptc_func:
                    ptc_func = _ptc_func[ptc.ptcFitType]
                    x = np.logspace(np.log10(np.nanmin(means)),
                                    np.log10(np.nanmax(means)), 100)
                    y = ptc_func(ptc.ptcFitPars[amp_name], x)
                    plt.plot(x, y, linestyle=':')
                plt.xscale('log')
                plt.yscale('log')
        plt.tight_layout(rect=(0.03, 0.03, 1, 0.95))
        fig.supxlabel('mean signal (ADU)')
        fig.supylabel('variance (ADU^2)')
        fig.suptitle(f'{det.getName()}')

        return pipeBase.Struct(ptc_plots=fig)


class PtcFpPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                dimensions=("instrument",)):
    ptc_results = cT.Input(name="ptc_results",
                           doc="PTC fit results",
                           storageClass="PhotonTransferCurveDataset",
                           dimensions=("instrument", "detector"),
                           isCalibration=True,
                           multiple=True,
                           deferLoad=True)
    ptc_a00 = cT.Output(name="ptc_a00_plot",
                        doc="Focal plane map of PTC a00",
                        storageClass="Plot",
                        dimensions=("instrument",))
    ptc_gain = cT.Output(name="ptc_gain_plot",
                         doc="Focal plane map of PTC gain",
                         storageClass="Plot",
                         dimensions=("instrument",))
    ptc_noise = cT.Output(name="ptc_noise_plot",
                          doc="Focal plane map of PTC noise",
                          storageClass="Plot",
                          dimensions=("instrument",))
    ptc_turnoff = cT.Output(name="ptc_turnoff_plot",
                            doc="Focal plane map of PTC turnoff",
                            storageClass="Plot",
                            dimensions=("instrument",))


class PtcFpPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                           pipelineConnections=PtcFpPlotsTaskConnections):
    """Configuration for PtcFpPlotsTask."""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)


class PtcFpPlotsTask(pipeBase.PipelineTask):
    """Create summary plots of the PTC results for the full focal plane."""
    ConfigClass = PtcFpPlotsTaskConfig
    _DefaultName = "ptcFpPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = (self.config.yfigsize, self.config.xfigsize)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        input_handles = butlerQC.get(inputRefs)

        ptc_results = {_.dataId['detector']: _ for
                       _ in input_handles['ptc_results']}

        struct = self.run(ptc_results)
        butlerQC.put(struct, outputRefs)
        plt.close()

    def run(self, ptc_results):
        """
        Create summary plots of PTC results for the full focal plane.

        Parameters
        ----------
        ptc_results : `dict`
            Dictionary of handles to `PhotonTransferCurveDataset`s,
            keyed by detector.

        Returns
        -------
        struct : `lsst.pipe.base.struct`
            Struct of `~matplotlib.figure.Figure`s keyed by PTC result name.
        """
        amp_data = defaultdict(lambda: defaultdict(dict))
        for detector, handle in ptc_results.items():
            ptc = handle.get()
            for amp, values in ptc.ptcFitPars.items():
                if len(values) != 3:
                    continue
                ptc_a00, ptc_gain, ptc_var = values
                amp_data['ptc_a00'][detector][amp] = -ptc_a00
                amp_data['ptc_gain'][detector][amp] = ptc_gain
                if ptc_var > 0:
                    amp_data['ptc_noise'][detector][amp] = np.sqrt(ptc_var)
                amp_data['ptc_turnoff'][detector][amp] \
                    = np.max(ptc.finalMeans[amp]) if ptc.finalMeans[amp] else -1

        plots = {}
        for field in amp_data:
            plots[field] = plt.figure(figsize=self.figsize)
            ax = plots[field].add_subplot(111)
            # TODO: figure out how to get the run number into the
            # plot title.
            plot_focal_plane(ax, amp_data[field], title=f"{field}",
                             z_range="clipped_autoscale")

        return pipeBase.Struct(**plots)


class RowMeansVarianceTaskConnections(pipeBase.PipelineTaskConnections,
                                      dimensions=("instrument", "detector")):

    ptc_results = cT.Input(
        name="ptc_results",
        doc="Photon Transfer Curve dataset",
        storageClass="PhotonTransferCurveDataset",
        dimensions=("instrument", "detector"),
        isCalibration=True)

    ptc_frames = cT.Input(
        name="ptc_frames",
        doc="The ISR'd flat pairs used in the PTC analysis",
        storageClass="Exposure",
        dimensions=("instrument", "exposure", "detector"),
        multiple=True,
        deferLoad=True)

    row_means_variance_plot = cT.Output(
        name="row_means_variance_plot",
        doc=("Plots of row_means_variance vs expected Poisson signal "
             "for each amp in a CCD for a PTC dataset"),
        storageClass="Plot",
        dimensions=("instrument", "detector"))

    row_means_variance_stats = cT.Output(
        name="row_means_variance_stats",
        doc=("Slopes of row_means_variance vs expected Poisson signal "
             "for each amp in a CCD for a PTC dataset"),
        storageClass="DataFrame",
        dimensions=("instrument", "detector"))


class RowMeansVarianceTaskConfig(pipeBase.PipelineTaskConfig,
                                 pipelineConnections=RowMeansVarianceTaskConnections):
    """Configuration for RowMeansVarianceTask"""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    nsig = pexConfig.Field(doc="Number of stdevs for variance clip",
                           dtype=float, default=10.0)
    min_signal = pexConfig.Field(doc="Minimum signal (e-/pixel) for fitting.",
                                 dtype=float, default=3000)
    max_signal = pexConfig.Field(doc="Maximum signal (e-/pixel) for fitting.",
                                 dtype=float, default=1e5)


class RowMeansVarianceTask(pipeBase.PipelineTask):
    """Task to compute row means variance for at PTC dataset."""
    ConfigClass = RowMeansVarianceTaskConfig
    _DefaultName = "rowMeansVarianceTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.nsig = self.config.nsig
        self.figsize = self.config.xfigsize, self.config.yfigsize
        self.min_signal = self.config.min_signal
        self.max_signal = self.config.max_signal

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        ptc_frames = {ref.dataId['exposure']: ref for
                      ref in inputs['ptc_frames']}
        # The expId pairs should be the same for all amps, so just get
        # them from one of the amps and strip out the extra layer of
        # Python list.
        ptc_results = inputs['ptc_results']
        expId_pairs = [_[0] for _ in
                       list(ptc_results.inputExpIdPairs.values())[0]]
        outputs = self.run(ptc_frames, expId_pairs, ptc_results.gain)
        butlerQC.put(outputs, outputRefs)

    def run(self, ptc_frames, expId_pairs, gains):
        data = defaultdict(list)
        for expIds in expId_pairs:
            frame1 = ptc_frames[expIds[0]].get()
            det = frame1.getDetector()
            flat1 = frame1.getMaskedImage()
            flat2 = ptc_frames[expIds[1]].get().getMaskedImage()
            det_name = det.getName()
            for amp in det:
                amp_name = amp.getName()
                bbox = amp.getBBox()
                stats1 = afwMath.makeStatistics(flat1[bbox], afwMath.MEAN)
                stats2 = afwMath.makeStatistics(flat2[bbox], afwMath.MEAN)
                signal = (stats1.getValue(afwMath.MEAN) +
                          stats2.getValue(afwMath.MEAN))/2.*gains[amp_name]
                diff = flat1.Factory(flat1, bbox, deep=True)
                diff -= flat2[bbox]
                row_means = np.mean(diff.getImage().array, axis=1,
                                    dtype=np.float64)
                sctrl = afwMath.StatisticsControl()
                sctrl.setNumSigmaClip(self.nsig)
                row_mean_variance = afwMath.makeStatistics(
                    row_means, afwMath.VARIANCECLIP, sctrl).getValue()
                row_mean_variance *= gains[amp_name]**2
                data['det_name'].append(det_name)
                data['amp_name'].append(amp_name)
                data['signal'].append(signal)
                data['row_mean_variance'].append(row_mean_variance)
        df0 = pd.DataFrame(data)

        fig = plt.figure(figsize=self.figsize)
        data = defaultdict(list)
        for amp in det:
            amp_name = amp.getName()
            numcols = amp.getBBox().width
            df = df0.query(f"amp_name == '{amp.getName()}'")
            signal = df['signal'].to_numpy()
            row_mean_var = df['row_mean_variance'].to_numpy()
            index = np.where((self.min_signal < signal)
                             & (signal < self.max_signal)
                             & (row_mean_var == row_mean_var))
            slope = sum(row_mean_var[index])/sum(2.*signal[index]/numcols)
            data['det_name'].append(det_name)
            data['amp_name'].append(amp_name)
            data['slope'].append(slope)
            plt.scatter(2*signal[index]/numcols, row_mean_var[index], s=2,
                        label=amp.getName())
        plt.xscale('log')
        plt.yscale('log')
        xmin, xmax, ymin, ymax = plt.axis()
        xymin = min(xmin, ymin)
        xymax = max(xmax, ymax)
        plt.plot([xymin, xymax], [xymin, xymax], linestyle=':')
        plt.legend(fontsize='x-small', ncol=2, loc=2)
        plt.axis((xmin, xmax, ymin, ymax))
        plt.xlabel('2*(flux/(e-/pixel))/num_cols')
        plt.ylabel('var(row_means)')
        plt.title(det_name)

        return pipeBase.Struct(row_means_variance_plot=fig,
                               row_means_variance_stats=pd.DataFrame(data))


class RowMeansVarianceFpPlotTaskConnections(pipeBase.PipelineTaskConnections,
                                            dimensions=("instrument",)):
    row_means_variance_stats = cT.Input(
        name="row_means_variance_stats",
        doc=("Slopes of row_means_variance vs expected Poisson signal "
             "for each amp in a CCD for a PTC dataset"),
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

    row_means_variance_slopes = cT.Output(
        name="row_means_variance_slopes_plot",
        doc="Focal plan map of slope of row means variance.",
        storageClass="Plot",
        dimensions=("instrument",))


class RowMeansVarianceFpPlotTaskConfig(pipeBase.PipelineTaskConfig,
                                       pipelineConnections=RowMeansVarianceFpPlotTaskConnections):
    """Configuration for RowMeansVarianceFpPlotTask."""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)


class RowMeansVarianceFpPlotTask(pipeBase.PipelineTask):
    """Create focal plane mosaic of slope of row means variance."""
    ConfigClass = RowMeansVarianceFpPlotTaskConfig
    _DefaultName = "rowMeansVarianceFpPlotTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.xfigsize, self.config.yfigsize

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        camera = inputs['camera']

        rmv_stats = {_.dataId['detector']: _ for _
                     in inputs['row_means_variance_stats']}

        struct = self.run(rmv_stats, camera)
        butlerQC.put(struct, outputRefs)
        plt.close()

    def run(self, rmv_stats, camera):
        amp_data = defaultdict(dict)
        for detector, ref in rmv_stats.items():
            df = ref.get()
            for _, row in df.iterrows():
                amp_data[row.det_name][row.amp_name] = row.slope
        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(111)
        plot_focal_plane(ax, amp_data, title="row means variance slopes",
                         camera=camera, z_range="clipped_autoscale")

        return pipeBase.Struct(row_means_variance_slopes=fig)
