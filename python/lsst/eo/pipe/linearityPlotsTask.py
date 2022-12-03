from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from lsst.eo.pipe.plotting import plot_focal_plane

__all__ = ['LinearityPlotsTask', 'LinearityFpPlotsTask']


def get_pd_values(pd_integrals, ptc, amp_name='C10'):
    values = []
    for pair in ptc.inputExpIdPairs[amp_name]:
        if pair[0][0] not in pd_integrals or pair[0][1] not in pd_integrals:
            # Use sentinel value of -1 to enable this entry to be deselected.
            values.append(-1)
        else:
            values.append((pd_integrals[pair[0][0]] +
                           pd_integrals[pair[0][1]])/2.)
    return np.array(values)


def linearity_fit(flux, Ne, y_range=(1e3, 9e4)):
    max_index = np.where(Ne == max(Ne))[0][0]
    index = np.where((y_range[0] < Ne) & (Ne < y_range[1])
                     & (flux <= flux[max_index]) & (flux > 0))
    aa = sum(flux[index]*Ne[index])/sum(flux[index]**2)

    def func(x):
        return aa*x

    resids = (Ne - func(flux))/Ne
    return func, resids, index


class LinearityPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                    dimensions=("instrument", "detector")):
    pd_data = cT.Input(
        name="photodiode",
        doc="Photodiode readings data.",
        storageClass="IsrCalib",
        dimensions=("instrument", "exposure"),
        multiple=True,
        deferLoad=True)

    ptc_results = cT.Input(
        name="ptc_results",
        doc="PTC results dataset",
        storageClass="PhotonTransferCurveDataset",
        dimensions=("instrument", "detector"),
        isCalibration=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    linearity_plots = cT.Output(
        name="linearity_fit_plot",
        doc="Plot of linearity fit to mean amp signal vs pd integral",
        storageClass="Plot",
        dimensions=("instrument", "detector"))

    residual_plots = cT.Output(
        name="linearity_residuals_plot",
        doc="Plot of linearity fit residuals vs pd integral",
        storageClass="Plot",
        dimensions=("instrument", "detector"))

    linearity_results = cT.Output(
        name="linearity_results",
        doc="Data frame of linearity results",
        storageClass="DataFrame",
        dimensions=("instrument", "detector"))


class LinearityPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                               pipelineConnections=LinearityPlotsTaskConnections):
    """Configuration for LinearityPlotsTask."""
    fit_range_min = pexConfig.Field(doc="Fit range lower bound in e-",
                                    dtype=float, default=1e3)
    fit_range_max = pexConfig.Field(doc="Fit range upper bound in e-",
                                    dtype=float, default=9e4)
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)


class LinearityPlotsTask(pipeBase.PipelineTask):
    """Create linearity plots from a flat pairs dataset using PTC results."""
    ConfigClass = LinearityPlotsTaskConfig
    _defaultName = "linearityPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.fit_range = self.config.fit_range_min, self.config.fit_range_max
        self.figsize = self.config.yfigsize, self.config.xfigsize

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        pd_data = {}
        for handle in inputs["pd_data"]:
            exposure = handle.dataId["exposure"]
            physical_filter = handle.dataId.records["exposure"].physical_filter
            pd_data[exposure] = (physical_filter, handle)
        acq_run = handle.dataId.records["exposure"].science_program

        ptc_ref = inputs['ptc_results']

        camera = inputs['camera']

        struct = self.run(pd_data, ptc_ref, camera, acq_run)
        butlerQC.put(struct, outputRefs)
        plt.close()

    def run(self, pd_data, ptc_ref, camera, acq_run):
        # Compute pd_integrals, applying Run 5 correction factors:
        def pd_corr(physical_filter):
            corr_factors = {'SDSSi': 1,
                            'SDSSi~ND_OD0.5': 0.9974873256470026}
            return corr_factors.get(physical_filter, 1)

        pd_integrals = {}
        for exposure, (physical_filter, handle) in pd_data.items():
            photodiode = handle.get()
            pd_integrals[exposure] \
                = photodiode.integrate()[0]*pd_corr(physical_filter)

        # Fit linear model to each amp.
        detector = ptc_ref.dataId['detector']
        ptc = ptc_ref.get()
        pd_values = get_pd_values(pd_integrals, ptc)
        det = camera[detector]
        det_name = det.getName()
        resids = {}
        index = {}
        linearity_plots = plt.figure(figsize=self.figsize)
        amp_data = defaultdict(list)
        for i, amp in enumerate(det, 1):
            amp_name = amp.getName()
            gain = amp.getGain()
            Ne = np.array(ptc.rawMeans[amp_name])*gain
            try:
                func, resids[amp_name], index[amp_name] \
                    = linearity_fit(pd_values, Ne, y_range=self.fit_range)
            except (IndexError, ZeroDivisionError) as eobj:
                pass
            else:
                # Results to include in final data frame:
                amp_data['det_name'].append(det_name)
                amp_data['amp_name'].append(amp_name)
                amp_data['max_frac_dev'].append(
                    max(np.abs(resids[amp_name][index[amp_name]])))
                amp_data['max_observed_signal'].append(
                    max(ptc.rawMeans[amp_name]))
                # plot for current amp
                ax = plt.subplot(4, 4, i)
                plt.scatter(pd_values, Ne, s=2)
                plt.scatter(pd_values[index[amp_name]], Ne[index[amp_name]],
                            s=2, color='red')
                plt.xscale('log')
                plt.yscale('log')
                plt.axhline(1e3, linestyle=':')
                plt.axhline(9e4, linestyle=':')
                xlims = ax.get_xlim()
                xvals = np.logspace(np.log10(xlims[0]), np.log10(xlims[1]), 100)
                plt.plot(xvals, func(xvals), linestyle='--')
        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.suptitle(f'linearity curves, run {acq_run}, {det_name}')

        residual_plots = plt.figure(figsize=self.figsize)
        for i, amp in enumerate(det, 1):
            ax = plt.subplot(4, 4, i)
            amp_name = amp.getName()
            try:
                plt.scatter(pd_values, resids[amp_name], s=2)
            except KeyError:
                pass
            else:
                plt.scatter(pd_values[index[amp_name]],
                            resids[amp_name][index[amp_name]],
                            s=2, color='red')
                plt.xscale('log')
                plt.axhline(0.02, linestyle=':')
                plt.axhline(-0.02, linestyle=':')
                plt.axhline(0, linestyle='--')
                plt.ylim(-0.03, 0.03)
        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.suptitle(f'linearity residuals, run {acq_run}, {det_name}')

        linearity_results = pd.DataFrame(amp_data)

        outputs = {'linearity_plots': linearity_plots,
                   'residual_plots': residual_plots,
                   'linearity_results': linearity_results}
        return pipeBase.Struct(**outputs)


class LinearityFpPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                      dimensions=("instrument",)):
    linearity_results = cT.Input(
        name="linearity_results",
        doc="Linearity results for each detector",
        storageClass="DataFrame",
        dimensions=("instrument", "detector"),
        multiple=True,
        deferLoad=True)

    max_frac_dev = cT.Output(name="max_frac_dev",
                             doc=("Maximum fractional deviation from "
                                  "linear fit of signal vs flux"),
                             storageClass="Plot",
                             dimensions=("instrument",))

    max_observed_signal = cT.Output(name="max_observed_signal",
                                    doc=("Maximum observed signal (ADU) "
                                         "from flat pair acquisition."),
                                    storageClass="Plot",
                                    dimensions=("instrument",))


class LinearityFpPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                                 pipelineConnections=LinearityFpPlotsTaskConnections):
    """Configuration for LinearityFpPlotsTask."""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)


class LinearityFpPlotsTask(pipeBase.PipelineTask):
    """Create summary plots of linearity analysis results for the
    full focal plane."""
    ConfigClass = LinearityFpPlotsTaskConfig
    _defaultName = "linearityFpPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = (self.config.yfigsize, self.config.xfigsize)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        handles = butlerQC.get(inputRefs)
        results = {_.dataId['detector']: _ for _
                   in handles["linearity_results"]}
        struct = self.run(results)
        butlerQC.put(struct, outputRefs)
        plt.close()

    def run(self, results):
        """
        Create summary plots of linearity analysis results for the full
        focal plane.

        Parameters
        ----------
        results : `dict`
            Dictionary of handles to `DataFrame`s containing the linearity
            results keyed by detector.

        Returns
        -------
        struct : `lsst.pipe.base.struct`
            Struct of `~matplotlib.figure.Figure`s keyed by result name.
        """
        amp_data = defaultdict(lambda: defaultdict(dict))
        for detector, handle in results.items():
            df = handle.get()
            columns = set(df.columns)
            columns.discard('amp_name')
            columns.discard('det_name')
            for _, row in df.iterrows():
                for column in columns:
                    amp_data[column][row.det_name][row.amp_name] = row[column]
        plots = {}
        for column in set(amp_data.keys()):
            plots[column] = plt.figure(figsize=self.figsize)
            ax = plots[column].add_subplot(111)
            plot_focal_plane(ax, amp_data[column], camera=camera,
                             z_range=None, title=f"{column}")

        return pipeBase.Struct(**plots)
