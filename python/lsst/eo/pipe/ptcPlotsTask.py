from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from lsst.eo.pipe.plotting import plot_focal_plane


__all__ = ['PtcPlotsTask', 'PtcFpPlotsTask']


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
            plt.scatter(means, variances, s=4, color='blue')
            plt.scatter(means[index], variances[index], s=2, color='green')
            if means.size != 0:
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
                               dtype=float, default=10)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=12)


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
