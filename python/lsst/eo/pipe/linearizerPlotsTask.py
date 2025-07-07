from collections import defaultdict
import numpy as np
import lsst.afw.image as afw_image
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
import matplotlib.pyplot as plt
from .plotting import plot_focal_plane, hist_amp_data, append_acq_run


__all__ = ['LinearizerFpPlotsTask']


def infer_max_frac_devs(linearizer, ptc, det, Ne_min=1e3, Ne_max=9e4):
    """
    Back out the maximum fractional deviation from a linear fit using
    the cp_pipe linearizer object for each amp in a detector.  This
    function assumes a spline model for the linearity correction was
    used.
    """
    max_frac_devs = {}
    for amp in det:
        amp_name = amp.getName()
        # Extract the measured ADU values and fit residuals.
        yvals = ptc.rawMeans[amp_name]
        # Replace nans in the fit residuals with zeros.
        resids = np.zeros(len(yvals))
        for i, value in enumerate(linearizer.fitResiduals[amp_name]):
            if value == value:
                resids[i] = value

        # Create an uncorrected image with a ramp covering the
        # measured ADU values from which to infer the measure signal
        # vs incident (corrected) flux (used as proxy for pd integrals).
        bbox = linearizer.linearityBBox[amp_name]
        image = afw_image.ImageF(bbox.width, bbox.height)
        # Subtract off the fit residuals so that measured points lie
        # on the spline-fitted linearity curve.
        meas_min = min(yvals - resids)
        meas_max = max(yvals - resids)
        meas = np.logspace(np.log10(meas_min)*np.ones(bbox.width),
                           np.log10(meas_max)*np.ones(bbox.width),
                           bbox.height,
                           axis=0)
        # Apply the linearity correction.
        image.array[:] = meas.copy()
        coeffs = linearizer.linearityCoeffs[amp_name]
        correction_func = linearizer.getLinearityTypeByName(
            linearizer.linearityType[amp_name])()
        gain = ptc.gain[amp_name]
        correction_func(image=image, coeffs=coeffs, gain=gain)
        corr = image.array.copy()
        # Evaluate the corrected flux values at the measured ADU values.
        xvals = np.interp(yvals - resids, meas[:, 0], corr[:, 0])
        # Select points within the desired fit range.
        index = np.where((1e3 < yvals*gain)*(yvals*gain < 9e4))
        # Fit line with zero y-intercept.
        denom = sum(xvals[index]**2/yvals[index])
        if denom == 0:
            continue
        aa = sum(xvals[index])/denom
        def func(x):
            return aa*x
        # Compute the fractional residuals to this fit.
        frac_resids = (yvals - func(xvals))/yvals
        max_frac_devs[amp_name] = float(max(np.abs(frac_resids)))

    return max_frac_devs

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
            max_frac_devs = infer_max_frac_devs(linearizer, ptc, det)
            for amp in det:
                amp_name = amp.getName()
                frac_defs = (linearizer.fitResiduals[amp_name]
                             / ptc.rawMeans[amp_name])
                try:
                    amp_data['max_frac_dev'][det_name][amp_name] \
                        = max_frac_devs[amp_name]
                except KeyError:
                    pass
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
