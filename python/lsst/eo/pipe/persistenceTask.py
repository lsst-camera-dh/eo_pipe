from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import lsst.afw.math as afw_math
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

from .plotting import append_acq_run
from .dsref_utils import get_plot_locations_by_dstype


__all__ = ['PersistenceTask']


def get_plot_locations(repo, collections):
    dstypes = ('persistence_plot',)
    return get_plot_locations_by_dstype(repo, collections, dstypes)


class PersistenceTaskConnections(pipeBase.PipelineTaskConnections,
                                 dimensions=("instrument", "detector")):
    exposure_handles = cT.Input(
        name="eoDarkIsr",
        doc="ISR'd dark frame to use for persistence measurement",
        storageClass="Exposure",
        dimensions=("instrument", "exposure", "detector"),
        multiple=True,
        deferLoad=True)

    persistence_data = cT.Output(
        name="persistence_data",
        doc=("Per-amp persistent signal measurements of dark frames following "
             "a high flux flat exposure."),
        storageClass="DataFrame",
        dimensions=("instrument", "detector"))

    persistence_plot = cT.Output(
        name="persistence_plot",
        doc=("Plot of per-amp persistent signal measurements as a function "
             "of time for dark frames following a high flux flat exposure."),
        storageClass="Plot",
        dimensions=("instrument", "detector"))


class PersistenceTaskConfig(pipeBase.PipelineTaskConfig,
                            pipelineConnections=PersistenceTaskConnections):
    nsig = pexConfig.Field(doc="Number of stdevs for variance clip",
                           dtype=float, default=5.0)
    niterclip = pexConfig.Field(doc="Number of iterations to use for "
                                "variance clip",
                                dtype=int, default=3)
    maskNameList = pexConfig.ListField(
        doc="Mask list to exclude from statistics calculations.",
        dtype=str,
        default=["SUSPECT", "BAD", "NO_DATA", "SAT"])
    acq_run = pexConfig.Field(doc="Acquisition run number.",
                              dtype=str, default="")


class PersistenceTask(pipeBase.PipelineTask):
    """
    Measure persistent signal in dark frames following a high flux
    flat exposure.
    """
    ConfigClass = PersistenceTaskConfig
    _DefaultName = "persistenceTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.nsig = self.config.nsig
        self.niter = self.config.niterclip
        self.masks = self.config.maskNameList

    def run(self, exposure_handles):
        data = defaultdict(list)
        flags = afw_math.MEANCLIP | afw_math.STDEVCLIP
        for handle in exposure_handles:
            exp = handle.get()
            md = exp.getMetadata()
            mjd_trg = md.get('MJD-TRG')
            image_type = handle.dataId.records['exposure'].observation_type
            mask_val = exp.getMask().getPlaneBitMask(self.masks)
            sctrl = afw_math.StatisticsControl(self.nsig, self.niter, mask_val)
            sctrl.setAndMask(mask_val)

            det = exp.getDetector()
            det_name = det.getName()
            for amp in det:
                amp_name = amp.getName()
                amp_image = exp.getImage()[amp.getBBox()]
                stats = afw_math.makeStatistics(amp_image, flags, sctrl)
                data['det_name'].append(det_name)
                data['amp_name'].append(amp_name)
                data['image_type'].append(image_type)
                data['tseqnum'].append(md.get('TSEQNUM'))
                data['mjd_trg'].append(mjd_trg)
                data['dark_time'].append(md.get('DARKTIME'))
                data['mean_signal'].append(stats.getValue(afw_math.MEANCLIP))
                data['stdev'].append(stats.getValue(afw_math.STDEVCLIP))
        df0 = pd.DataFrame(data)
        mjd0 = min(data['mjd_trg'])
        df0['time'] = (df0['mjd_trg'] - mjd0)*24.0*60.0

        fig = plt.figure()
        for amp in det:
            amp_name = amp.getName()
            df = df0.query(f"amp_name=='{amp_name}' and image_type=='dark'")
            plt.scatter(df['time'], df['mean_signal'], s=2, label=amp_name)
        plt.legend(fontsize='x-small', ncol=2)
        # Include 0 in x-axis range.
        axis = plt.axis()
        plt.xlim(0, axis[1])
        plt.xlabel('Time since flat readout trigger (min)')
        plt.ylabel('Mean residual signal (ADU)')
        plt.title(append_acq_run(self, "Persistence test", det_name))
        return pipeBase.Struct(persistence_data=df0,
                               persistence_plot=fig)
