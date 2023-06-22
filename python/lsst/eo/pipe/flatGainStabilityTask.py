from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import lsst.afw.math as afw_math
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from .isr_utils import apply_minimal_isr
from .dsref_utils import get_plot_locations_by_dstype
from .plotting import append_acq_run


__all__ = ['FlatGainStabilityTask', 'FlatGainStabilityPlotsTask']


def get_plot_locations(repo, collections):
    dstypes = ('flat_gain_stability_plots',)
    return get_plot_locations_by_dstype(repo, collections, dstypes)


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

    pd_data = cT.Input(
        name="photodiode",
        doc="Photodiode readings data.",
        storageClass="IsrCalib",
        dimensions=("instrument", "exposure"),
        multiple=True,
        deferLoad=True)

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
    pd_integration_method = pexConfig.ChoiceField(
        dtype=str,
        doc="Integration method for photodiode monitoring data.",
        default="CHARGE_SUM",
        allowed={
            "DIRECT_SUM": ("Use numpy's trapz integrator on all photodiode "
                           "readout entries"),
            "TRIMMED_SUM": ("Use numpy's trapz integrator, clipping the "
                            "leading and trailing entries, which are "
                            "nominally at zero baseline level."),
            "CHARGE_SUM": ("Treat the current values as integrated charge "
                           "over the sampling interval and simply sum "
                           "the values, after subtracting a baseline level."),
        })
    pd_current_scale = pexConfig.Field(
        dtype=float,
        doc="Scale factor to apply to photodiode current values for the "
            "``CHARGE_SUM`` integration method.",
        default=-1.0,
    )
    ccob_led_constraint = pexConfig.Field(
        dtype=str,
        doc="Value of CCOBLED keyword to use for filtering sflat exposures. "
        "If 'None', then process all of the exposures.",
        default="None")


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
        self.pd_integration_method = self.config.pd_integration_method
        self.pd_current_scale = self.config.pd_current_scale
        self.ccob_led_constraint = self.config.ccob_led_constraint

    def run(self, raws, bias, dark, pd_data, camera):
        # Prepare dict of photodiode integrals and key by exposure.
        pd_integrals = {}
        for ref in pd_data:
            exposure = ref.dataId['exposure']
            pd_calib = ref.get()
            pd_calib.integrationMethod = self.pd_integration_method
            pd_calib.currentScale = self.pd_current_scale
            pd_integrals[exposure] = pd_calib.integrate()

        det = camera[raws[0].dataId['detector']]
        det_name = det.getName()
        raft, slot = det_name.split('_')
        data = defaultdict(list)
        num_raws = len(raws)
        self.log.info(f"Analyzing {det_name}, with {num_raws} flats")
        for i, handle in enumerate(raws):
            raw = handle.get()
            md = raw.getMetadata()
            exposure = handle.dataId['exposure']
            if (self.ccob_led_constraint != "None"
                and md.get('CCOBLED') != self.ccob_led_constraint):
                self.log.info(f"Skipping exposure {exposure} "
                              "because of ccob_led_constraint.")
                continue
            self.log.info(f"Processing exposure {exposure}, "
                          f"image {i} of {num_raws}")
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
                data['pd_integral'].append(pd_integrals[exposure])
                data['mjd'].append(mjd_obs)
                data['amp_name'].append(amp_name)
                data['median'].append(stats.getValue(afw_math.MEDIAN))
                data['stdev_clip'].append(stats.getValue(afw_math.STDEVCLIP))
        return pipeBase.Struct(flat_gain_stability_stats=pd.DataFrame(data))


class FlatGainStabilityPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                            dimensions=("instrument",)):
    flat_gain_stability_stats = cT.Input(
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
                               dtype=float, default=18)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=18)
    y_range_min = pexConfig.Field(doc="Minimum of y-axis range in each plot",
                                  dtype=float, default=0.998)
    y_range_max = pexConfig.Field(doc="Maximum of y-axis range in each plot",
                                  dtype=float, default=1.002)
    ccob_led_constraint = pexConfig.Field(
        dtype=str,
        doc="Value of CCOBLED keyword to use for filtering sflat exposures. "
        "If 'None', then process all of the exposures.",
        default="None")
    acq_run = pexConfig.Field(doc="Acquisition run number.",
                              dtype=str, default="")


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
        self.y_range = self.config.y_range_min, self.config.y_range_max
        self.ccob_led = self.config.ccob_led_constraint

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
                signal = stats['median'].to_numpy()
                amp_names = stats['amp_name'].to_numpy()
                pd_integrals = stats['pd_integral'].to_numpy()
                for amp in det:
                    index = np.where(amp_names == amp.getName())
                    pd_values = pd_integrals[index].copy()
                    pd_values /= np.mean(pd_values)
                    signal[index] /= pd_values
                    signal[index] /= np.mean(signal[index])
                plt.scatter(time, signal, s=2, label=slot)
            plt.ylim(*self.y_range)
            plt.legend(fontsize='x-small', ncol=2)
            plt.title(raft, fontsize='small')
        plt.tight_layout(rect=(0.03, 0.03, 1, 0.95))
        fig.supxlabel(f'(MJD - {int(t0)})*24 (hours)')
        fig.supylabel('counts/mean(counts)')
        if self.ccob_led != "None":
            suptitle = append_acq_run(self, 'flat gain stability',
                                      f"{self.ccob_led} LED")
        else:
            suptitle = append_acq_run(self, 'flat gain stability')
        fig.suptitle(suptitle)
        return pipeBase.Struct(flat_gain_stability_plots=fig)
