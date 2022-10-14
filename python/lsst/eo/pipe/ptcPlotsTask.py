from collections import defaultdict

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from lsst.eo.pipe.plotting import plot_focal_plane

class PtcPlotsTaskConnections(pipeBase.PipelineTaskCponnections,
                              dimensions=("instrument",)):
    ptc_results = cT.Input(name="ptc_results",
                           doc="PTC fit results",
                           storageClass="PhotonTransferCurveDataset",
                           dimensions=("instrument", "detector"),
                           multiple=True,
                           deferLoad=True)
    ptc_plots = cT.Output(name="ptc_plots",
                          doc="PNG plots of PTC fit results",
                          storageClass="Plot",
                          dimensions=("instrument",))

class PtcPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                         pipelineConnections=PtcPlotsTaskConnections):
    figsize = pexConfig.Field(doc="Figure size in inches.",
                              dtype=tuple, default=(12, 10))

class PtcPlotsTask(pipeBase.PipelineTask):

    ConfigClass = PtcPlotsTaskConfig
    _defaultName = "ptcPlotsTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.figsize

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        input_handles = butlerQC.get(inputRefs)

        ptc_results = {_.dataId['detector']: _ for
                       _ in input_handles['ptc_results']}

        struct = self.run(ptc_results)
        butlerQC.put(struct.ptc_plots, outputRefs.ptc_plots)

    def run(self, ptc_results):
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

