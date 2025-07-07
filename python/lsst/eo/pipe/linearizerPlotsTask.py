from collections import defaultdict
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
import matplotlib.pyplot as plt
from .plotting import plot_focal_plane, hist_amp_data, append_acq_run


__all__ = ['LinearizerFpPlotsTask']


class LinearizerFpPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                       dimensions=("instrument",)):
    linearizers = cT.Input(
        name="linearizers",
        doc="linearizer fit results datasets",
        storageClass="Linearizer",
        dimensions=("instrument", "detector"),
        isCalibration=True,
        multiple=True,
        deferLoad=True)

    linearizerPtcs = cT.Input(
        name="linearizerPtcs",
        doc="Input linearizer PTC datasets",
        storageClass="PhotonTransferCurveDataset",
        dimensions=("instrument", "detector"),
        isCalibration=True,
        multiple=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",))

    max_frac_dev = cT.Output(
        name="max_frac_dev",
        doc=("Maximum fractional deviation from "
             "linear fit of signal vs flux "
             "computed by cpLinearizer task"),
        storageClass="Plot",
        dimensions=("instrument",))

    max_frac_dev_hist = cT.Output(
        name="max_frac_dev_hist",
        doc=("Maximum fractional deviation from "
             "linear fit of signal vs flux "
             "computed by cpLinearizer task"),
        storageClass="Plot",
        dimensions=("instrument",))

    max_observed_signal = cT.Output(
        name="max_observed_signal",
        doc=("Maximum observed signal computed by cpLinearizer task"),
        storageClass="Plot",
        dimensions=("instrument",))

    max_observed_signal_hist = cT.Output(
        name="max_observed_signal_hist",
        doc=("Maximum observed signal computed by cpLinearizer task"),
        storageClass="Plot",
        dimensions=("instrument",))

    linearity_turnoff = cT.Output(
        name="linearity_turnoff",
        doc=("Linearity turnoff computed by cpLinearizer task"),
        storageClass="Plot",
        dimensions=("instrument",))

    linearity_turnoff_hist = cT.Output(
        name="linearity_turnoff_hist",
        doc=("Linearity turnoff computed by cpLinearizer task"),
        storageClass="Plot",
        dimensions=("instrument",))


class LinearizerFpPlotsTaskConfig(
        pipeBase.PipelineTaskConfig,
        pipelineConnections=LinearizerFpPlotsTaskConnections):
    """Configuration for LinearityFpPlotsTask."""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    acq_run = pexConfig.Field(doc="Acquisition run number.",
                              dtype=str, default="")


class LinearizerFpPlotsTask(pipeBase.PipelineTask):
    ConfigClass = LinearizerFpPlotsTaskConfig
    _DefaultName = "linearizerFpPlotsTask"

    _z_range = {'max_frac_dev': (0, 0.05),
                'max_observed_signal': (5e4, 2e5),
                'linearity_turnoff': (5e4, 2e5)}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = (self.config.yfigsize, self.config.xfigsize)

    def run(self, linearizers, linearizerPtcs, camera):
        # Extract the per-amp data from the cpLinearizer task outputs.
        linearizer_map = {_.dataId['detector']: _ for _ in linearizers}
        ptc_map = {_.dataId['detector']: _ for _ in linearizerPtcs}
        amp_data = defaultdict(lambda: defaultdict(dict))
        for detector in linearizer_map:
            det = camera[detector]
            det_name = det.getName()
            linearizer = linearizer_map[detector].get()
            ptc = ptc_map[detector].get()
            for amp in det:
                amp_name = amp.getName()
                frac_defs = (linearizer.fitResiduals[amp_name]
                             / ptc.rawMeans[amp_name])
                amp_data['max_frac_dev'][det_name][amp_name] \
                    = float(max(frac_defs))
                amp_data['max_observed_signal'][det_name][amp_name] \
                    = float(linearizer.linearityMaxSignal[amp_name])
                amp_data['linearity_turnoff'][det_name][amp_name] \
                    = float(linearizer.linearityTurnoff[amp_name])

        # Create the focal plane plots and histograms.
        plots = {}
        hists = {}
        for column in set(amp_data.keys()):
            plots[column] = plt.figure(figsize=self.figsize)
            ax = plots[column].add_subplot(111)
            title = append_acq_run(self, column)
            z_range = self._z_range.get(column, None)
            plot_focal_plane(ax, amp_data[column], camera=camera,
                             z_range=z_range, title=title)
            hists[f"{column}_hist"] = plt.figure()
            hist_amp_data(amp_data[column], column, hist_range=z_range,
                          title=title)
        return pipeBase.Struct(**plots, **hists)
