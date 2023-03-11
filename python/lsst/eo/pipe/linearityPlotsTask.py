from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.daf.butler as daf_butler
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from .plotting import plot_focal_plane, append_acq_run
from .dsref_utils import get_plot_locations_by_dstype


__all__ = ['LinearityPlotsTask', 'LinearityFpPlotsTask']


def get_amp_data(repo, collections):
    """Get linearity results for each amp in the focal plane."""
    butler = daf_butler.Butler(repo, collections=collections)
    dsrefs = list(set(butler.registry.queryDatasets('linearity_results',
                                                    findFirst=True)))
    amp_data = defaultdict(lambda: defaultdict(dict))
    fields = 'max_frac_dev max_observed_signal linearity_turnoff'.split()
    for dsref in dsrefs:
        df = butler.getDirect(dsref)
        for _, row in df.iterrows():
            for field in fields:
                amp_data[field][row.det_name][row.amp_name] = row[field]
    return {field: dict(data) for field, data in amp_data.items()}


def get_plot_locations(repo, collections):
    dstypes = ('linearity_fit_plot', 'linearity_residuals_plot',
               'max_frac_dev', 'max_observed_signal', 'linearity_turnoff')
    return get_plot_locations_by_dstype(repo, collections, dstypes)


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


def linearity_fit(flux, Ne, y_range=(1e3, 9e4), max_frac_dev=0.05, logger=None):
    """
    Fit a line with the y-intercept fixed to zero, using the
    signal counts Ne as the variance in the chi-square, i.e.,

    chi2 = sum( (Ne - aa*flux)**2/Ne )

    Minimizing chi2 wrt aa, gives

    aa = sum(flux) / sum(flux**2/Ne)

    Also apply the y_range selection to the signal counts, Ne,
    and omit any flux values above the Ne peak and any non-positive
    flux values.
    """
    max_index = np.where(Ne == max(Ne))[0][0]
    index = np.where((y_range[0] < Ne) & (Ne < y_range[1])
                     & (flux <= flux[max_index]) & (flux > 0))
    aa = sum(flux[index])/sum(flux[index]**2/Ne[index])

    def func(x):
        return aa*x

    resids = (Ne - func(flux))/Ne

    # Compute maximum signal consistent with the linear fit within
    # the specified maximum fractional deviation.
    try:
        linearity_turnoff = np.max(Ne[np.where(np.abs(resids) < max_frac_dev)])
    except ValueError:
        if logger is not None:
            logger.info("ValueError from linearity_turnoff calculation. "
                        "Setting to zero")
        linearity_turnoff = 0.0

    return func, resids, index, linearity_turnoff


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
    max_frac_dev_spec = pexConfig.Field(doc=("Maximum fractional deviation "
                                             "specification to use for "
                                             "linearity turnoff calculation."),
                                        dtype=float, default=0.05)
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)


class LinearityPlotsTask(pipeBase.PipelineTask):
    """Create linearity plots from a flat pairs dataset using PTC results."""
    ConfigClass = LinearityPlotsTaskConfig
    _DefaultName = "linearityPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.fit_range = self.config.fit_range_min, self.config.fit_range_max
        self.max_frac_dev_spec = self.config.max_frac_dev_spec
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
            gain = ptc.gain[amp_name]
            Ne = np.array(ptc.rawMeans[amp_name])*gain
            try:
                func, resids[amp_name], index[amp_name], linearity_turnoff \
                    = linearity_fit(pd_values, Ne, y_range=self.fit_range,
                                    max_frac_dev=self.max_frac_dev_spec,
                                    logger=self.log)
                linearity_turnoff /= gain  # Convert from e- to ADU
            except (IndexError, ZeroDivisionError):
                pass
            else:
                # Results to include in final data frame:
                amp_data['det_name'].append(det_name)
                amp_data['amp_name'].append(amp_name)
                amp_data['max_frac_dev'].append(
                    max(np.abs(resids[amp_name][index[amp_name]])))
                amp_data['max_observed_signal'].append(
                    max(ptc.rawMeans[amp_name]))
                amp_data['linearity_turnoff'].append(linearity_turnoff)
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
        plt.tight_layout(rect=(0.03, 0.03, 1, 0.95))
        linearity_plots.supxlabel('photodiode current integral')
        linearity_plots.supylabel('e-/pixel')
        linearity_plots.suptitle(f'Linearity Curves, acq. run {acq_run}, '
                                 f'{det_name}')

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
        plt.tight_layout(rect=(0.03, 0.03, 1, 0.95))
        residual_plots.supxlabel('photodiode current integral')
        residual_plots.supylabel('e-/pixel')
        residual_plots.suptitle(f'Linearity Residuals, acq. run {acq_run}, '
                                f'{det_name}')

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

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

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

    linearity_turnoff = cT.Output(name="linearity_turnoff",
                                  doc=("Maximum signal (ADU) consistent "
                                       "with the linearity fit within "
                                       "the maximum fractional deviation "
                                       "spec, nominally 0.05."),
                                  storageClass="Plot",
                                  dimensions=("instrument",))


class LinearityFpPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                                 pipelineConnections=LinearityFpPlotsTaskConnections):
    """Configuration for LinearityFpPlotsTask."""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    acq_run = pexConfig.Field(doc="Acquistion run number.",
                              dtype=str, default="")


class LinearityFpPlotsTask(pipeBase.PipelineTask):
    """Create summary plots of linearity analysis results for the
    full focal plane."""
    ConfigClass = LinearityFpPlotsTaskConfig
    _DefaultName = "linearityFpPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = (self.config.yfigsize, self.config.xfigsize)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        results = {_.dataId['detector']: _ for _
                   in inputs["linearity_results"]}
        camera = inputs['camera']
        struct = self.run(results, camera)
        butlerQC.put(struct, outputRefs)
        plt.close()

    def run(self, results, camera):
        """
        Create summary plots of linearity analysis results for the full
        focal plane.

        Parameters
        ----------
        results : `dict`
            Dictionary of handles to `DataFrame`s containing the linearity
            results keyed by detector.
        camera: `lsst.afw.cameraGeom.Camera`
            Camera object.

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
                             z_range=None, title=append_acq_run(self, column))

        return pipeBase.Struct(**plots)
