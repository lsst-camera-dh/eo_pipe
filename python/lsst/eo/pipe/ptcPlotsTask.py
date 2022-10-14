from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from lsst.eo.pipe.plotting import plot_focal_plane

__all__ = ['PtcPlotsTask']

class PtcPlotsTaskConnections(pipeBase.PipelineTaskConnections,
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

class PtcPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                         pipelineConnections=PtcPlotsTaskConnections):
    """Configuration for PtsPlotsTask."""
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=10)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=12)

class PtcPlotsTask(pipeBase.PipelineTask):
    """Create summary plots of the PTC results for the full focal plane."""
    ConfigClass = PtcPlotsTaskConfig
    _defaultName = "ptcPlotsTask"

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
        amp_data = defaultdict(lambda : defaultdict(dict))
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
            plot_focal_plane(ax, amp_data[field], title=f"{field}")

        return pipeBase.Struct(**plots)
