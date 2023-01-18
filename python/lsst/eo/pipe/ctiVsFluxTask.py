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
from .eperTask import compute_ctis


__all__ = ['CtiVsFluxTask']


def find_flat_pairs(raws, flux_keyword=""):
    # Make sure raws can be indexed.
    raws = list(raws)

    # Test if flux_keyword exists by checking first raw exposure.  If not,
    # then set to flux_keyword=None and use exposure time instead.
    if flux_keyword not in raws[0].get().getMetadata():
        flux_keyword = ""

    # Assemble the list of key values for each raw.
    key_values = []
    for handle in raws:
        if flux_keyword != "":
            key_values.append(handle.get().getMetadata()[flux_keyword])
        else:
            key_values.append(handle.dataId.records['exposure'].exposure_time)

    # Find the pairs based on the key_values.
    raw_pairs = []
    for unique_key in set(key_values):
        i0 = key_values.index(unique_key)
        try:
            i1 = key_values.index(unique_key, i0)
        except ValueError:
            # A complete key pair doesn't exist, so omit the raw exposure
            # associated with this key_value.
            pass
        else:
            raw_pairs.append((raws[i0], raws[i1]))

    return raw_pairs


class CtiVsFluxTaskConnections(pipeBase.PipelineTaskConnections,
                               dimensions=("instrument", "detector")):
    raws = cT.Input(
        name="raw",
        doc="Raw pixel data from flat pair dataset.",
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

    cti_vs_flux = cT.Output(
        name="cti_vs_flux",
        doc="Serial and parallel CTI values as a function of flux.",
        storageClass="DataFrame",
        dimensions=("instrument", "detector"))

    scti_vs_flux_plot = cT.Output(
        name="scti_vs_flux_plot",
        doc="Plot of serial CTI vs flux.",
        storageClass="Plot",
        dimensions=("instrument", "detector"))

    pcti_vs_flux_plot = cT.Output(
        name="scti_vs_flux_plot",
        doc="Plot of parallel CTI vs flux.",
        storageClass="Plot",
        dimensions=("instrument", "detector"))


class CtiVsFluxTaskConfig(pipeBase.PipelineTaskConfig,
                          pipelineConnections=CtiVsFluxTaskConnections):
    nx_skip = pexConfig.Field(
        doc=("Number columns at the leading and trailing edges of "
             "the serial overscan to omit when estimating the "
             "serial overscan correction."),
        default=4,
        dtype=int)
    overscan_pixels = pexConfig.Field(
        doc=("Number of overscan rows or columns to use for "
             "evaluating the trailed signal in the overscan regions."),
        default=3,
        dtype=int)
    oscan_method = pexConfig.ChoiceField(
        doc="Overscan modeling method",
        default="median_per_row",
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
    flux_keyword = pexConfig.Field(
        doc=("FITS header keyword with target flux value "
             "for finding image pairs"),
        default="",
        dtype=str)


class CtiVsFluxTask(pipeBase.PipelineTask):
    """Task to measure serial and parallel CTI as a function of incident flux
    using the EPER method."""
    ConfigClass = CtiVsFluxTaskConfig
    _DefaultName = "ctiVsFluxTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.dx = self.config.nx_skip
        self.npix = self.config.overscan_pixels
        self.oscan_method = self.config.oscan_method
        self.deg = self.config.polynomial_degree
        self.flux_keyword = self.config.flux_keyword

    def run(self, raws, bias, dark, camera):
        det = camera[raws[0].dataId['detector']]
        det_name = det.getName()

        # Compute the serial and parallel CTIs over all flat pairs and amps.
        data = defaultdict(list)
        flat_pairs = find_flat_pairs(raws, flux_keyword=self.flux_keyword)
        for flat0, flat1 in flat_pairs:
            raw0 = flat0.get()
            raw1 = flat1.get()
            for amp, amp_info in enumerate(det):
                amp_name = amp_info.getName()
                images = [apply_minimal_isr(raw, bias, dark, amp, dx=self.dx,
                                            oscan_method=self.oscan_method,
                                            deg=self.deg)
                          for raw in (raw0, raw1)]
                median_flat = afw_math.statisticsStack(images, afw_math.MEDIAN)
                signal = np.median(median_flat.array)
                scti, pcti = compute_ctis(median_flat, det[amp],
                                          npix=self.npix)
                data['det_name'].append(det_name)
                data['amp_name'].append(amp_name)
                data['signal'].append(signal)
                data['scti'].append(scti)
                data['pcti'].append(pcti)
        df0 = pd.DataFrame(data)

        # Serial CTI vs flux plot.
        scti_vs_flux_plot = plt.figure()
        for amp in det:
            amp_name = amp.getName()
            df = df0.query(f"amp_name=='{amp_name}'")
            plt.scatter(df['signal'], df['scti'], label=amp_name)
        plt.legend(fontsize='x-small', ncol=4)
        plt.xlabel("Signal [ADU]")
        plt.ylabel("Serial CTI")
        plt.xscale("log")
        plt.yscale("log")
        plt.title(f"Serial CTI, {det_name}")

        # Parallel CTI vs flux plot.
        pcti_vs_flux_plot = plt.figure()
        for amp in det:
            amp_name = amp.getName()
            df = df0.query(f"amp_name=='{amp_name}'")
            plt.scatter(df['signal'], df['pcti'], label=amp_name)
        plt.legend(fontsize='x-small', ncol=4)
        plt.xlabel("Signal [ADU]")
        plt.ylabel("Parallel CTI")
        plt.xscale("log")
        plt.yscale("log")
        plt.title(f"Parallel CTI, {det_name}")

        return pipeBase.Struct(cti_vs_flux=df0,
                               scti_vs_flux_plot=scti_vs_flux_plot,
                               pcti_vs_flux_plot=pcti_vs_flux_plot)
