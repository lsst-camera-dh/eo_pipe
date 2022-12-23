from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import lsst.afw.math as afw_math
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from .isr_utils import apply_minimal_isr


class FlatGainStabilityTaskConnections(pipeBase.PipelineTaskConnections,
                                       dimensions=("instrument", "detector")):
    raws = cT.Input(
        name="raw",
        doc="Raw pixel data taken for flat gain stability tests.",
        storageClass="Exposure",
        dimensions=("instrument", "detector", "exposure"),
        multiple=True,
        deferLoad=True)

    bias = cT.Input(
        name="bias_frame",
        doc="Combined bias frame",
        storageClass="Exposure",
        dimensions=("instrument", "detector"),
        isCalibration=True)

    dark = cT.Input(
        name="dark_frame",
        doc="Combined dark frame",
        storageClass="Exposure",
        dimensions=("instrument", "detector"),
        isCalibration=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    flat_gain_stability_stats = cT.Output(
        name="flat_gain_stability_stats",
        doc=("Flat median signal for each amp, normalized by photodiode "
             "integrals, as a function of time."),
        storageClass="DataFrame",
        dimensions=("instrument", "detector"))


class FlatGainStabilityTaskConfig(pipeBase.PipelineTaskConfig,
                                  pipelineConnections=FlatGainStabilityTaskConnections):
    nx_skip = pexConfig.Field(
        doc=("Number columns at the leading and trailing edges of "
             "the serial overscan to omit when estimating the "
             "serial overscan correction."),
        default=4,
        dtype=int)
    oscan_method = pexConfig.ChoiceField(
        doc="Overscan modeling method",
        default="1d_poly",
        dtype=str,
        allowed={
            "mean": "Mean of all selected pixels in overscan region",
            "median": "Median of all selected pixels in overscan region",
            "median_per_row": "Median of each row of selected pixels",
            "1d_poly": "1D polynomial of degree 2 fit to median_per_row data"})
    polynomial_degree = pexConfig.Field(
        doc="Degree of polynomial to fit to overscan row medians",
        default=2,
        dtype=int)


class FlatGainStabilityTask(pipeBase.PipelineTask):
    """
    Task to measure flat gain stability statistics for a series of
    exposures of a CCD.
    """
    ConfigClass = FlatGainStabilityTaskConfig
    _DefaultName = "flatGainStabilityTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.dx = self.config.nx_skip
        self.oscan_method = self.config.oscan_method
        self.deg = self.config.polynomial_degree

    def run(self, raws, bias, dark, camera):
        det = camera[raws[0].dataId['detector']]
        det_name = det.getName()
        raft, slot = det_name.split('_')
        data = defaultdict(list)
        for handle in raws:
            raw = handle.get()
            md = raw.getMetadata()
            exposure = handle.dataId['exposure']
            mjd_obs = md.get('MJD-OBS')
            for amp, amp_info in enumerate(det):
                amp_name = amp_info.getName()
                image = apply_minimal_isr(raw, bias, dark, amp, dx=self.dx,
                                          oscan_method=self.oscan_method,
                                          deg=self.deg)
                flags = afw_math.MEDIAN | afw_math.STDEVCLIP
                stats = afw_math.makeStatistics(image, flags)
                data['raft'].append(raft)
                data['slot'].append(slot)
                data['exposure'].append(exposure)
                data['mjd'].append(mjd_obs)
                data['amp_name'].append(amp_name)
                data['median'].append(stats.getValue(afw_math.MEDIAN))
                data['stdev_clip'].append(stats.getValue(afw_math.STDEVCLIP))
        return pipeBase.Struct(flat_gain_stability_stats=pd.DataFrame(data))


class FlatGainStabilityPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                            dimensions=("instrument")):
    flat_gain_stabilty_stats = cT.Input(
        name="flat_gain_stability_stats",
        doc=("Flat median signal for each amp, normalized by photodiode "
             "integrals, as a function of time."),
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

    flat_gain_stability_plots = cT.Output(
        name="flat_gain_stability_plots",
        doc="Plots of median signal vs time for flat gain stability data",
        storageClass="Plot",
        dimensions=("instrument",))


class FlatGainStabilityPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                                       pipelineConnections=FlatGainStabilityPlotsTaskConnections):
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)


class FlatGainStabilityPlotsTask(pipeBase.PipelineTask):
    """
    Task to plot flat gain stability as a function of time for
    all of the detectors in the focal plane.
    """
    ConfigClass = FlatGainStabilityPlotsTaskConfig
    _DefaultName = "flatGainStabilityPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.yfigsize, self.config.xfigsize

    def run(self, flat_gain_stability_stats, camera):
        # Sort by raft and slot:
        handles = defaultdict(dict)
        for ref in flat_gain_stability_stats:
            det = camera[ref.dataId['detector']]
            raft, slot = det.getName().split('_')
            handles[raft][slot] = ref
        fig = plt.figure(figsize=self.figsize)
        rafts = sorted(handles.keys())
        t0 = None
        for i, raft in enumerate(rafts, 1):
            plt.subplot(5, 5, i)
            slots = sorted(handles[raft].keys())
            for slot in slots:
                stats = handles[raft][slot].get()
                if t0 is None:
                    t0 = np.floor(min(stats['mjd']))
                time = 24.*(stats['mjd'].to_numpy() - t0)
                signal = stats['median']/np.mean(stats['median'])
                plt.scatter(time, signal, s=2, label=slot)
            plt.legend(fontsize='x-small', ncol=2)
            plt.title(raft)
        plt.tight_layout(rect=(0.03, 0.03, 1, 0.95))
        fig.supxlabel(f'(MJD - {int(t0)}*24 (hours)')
        fig.supylabel('counts/mean(counts)')
        return pipeBase.Struct(flat_gain_stability_plots=fig)
