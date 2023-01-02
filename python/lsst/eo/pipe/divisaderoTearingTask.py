from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.stats
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.daf.butler as daf_butler
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from lsst.eo.pipe.plotting import plot_focal_plane


__all__ = ['DivisaderoTearingTask', 'DivisaderoFpPlotsTask']


def get_amp_data(repo, collections):
    """Get divisadero tearing results for each amp in the focal plane."""
    butler = daf_butler.Butler(repo, collections=collections)
    dsrefs = list(set(butler.registry.queryDatasets('divisadero_stats',
                                                    findFirst=True)))
    field = 'divisadero_tearing'
    amp_data = defaultdict(dict)
    for dsref in dsrefs:
        df = butler.getDirect(dsref)
        for _, row in df.iterrows():
            amp_data[row.det_name][row.amp_name] = row[field]
    return {field: dict(amp_data)}


class DivisaderoTearingTaskConnections(pipeBase.PipelineTaskConnections,
                                       dimensions=("instrument", "detector")):
    flat = cT.Input(
        name="flat",
        doc="Input combined flat to process",
        storageClass="Exposure",
        dimensions=("band", "instrument", "detector", "physical_filter"),
        isCalibration=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    divisadero_response = cT.Output(
        name="divisadero_response",
        doc="mean divisadero response vs column for each CCD half",
        storageClass="NumpyArray",
        dimensions=("instrument", "detector"))

    divisadero_stats = cT.Output(
        name="divisadero_stats",
        doc=("Divisadero statistics for all CCD amp boundaries "
             "in the focal plane"),
        storageClass="DataFrame",
        dimensions=("instrument", "detector"))


class DivisaderoTearingTaskConfig(pipeBase.PipelineTaskConfig,
                                  pipelineConnections=DivisaderoTearingTaskConnections):
    """Configuration for divisadero tearing task."""
    row_min = pexConfig.Field(doc="Minimum row to consider.",
                              dtype=int, default=10)
    row_max = pexConfig.Field(doc="Maximum row to consider.",
                              dtype=int, default=210)
    nedge = pexConfig.Field(doc=("Number of pixels to omit at edges of "
                                 "segments for removing the linear trend from "
                                 "the divisadero profiles vs column."),
                            dtype=int, default=25)
    npix = pexConfig.Field(doc=("Number of pixels inset from segment edges "
                                "to consider for divisadero tearing."),
                           dtype=int, default=2)


class DivisaderoTearingTask(pipeBase.PipelineTask):
    """
    Run the divisadero tearing analysis on combined flats over all of
    the CCDs in the focal plane.
    """
    ConfigClass = DivisaderoTearingTaskConfig
    _DefaultName = 'divisaderoTearingTask'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.row_min = self.config.row_min
        self.row_max = self.config.row_max
        self.nedge = self.config.nedge
        self.npix = self.config.npix

    def run(self, flat, camera):
        response_vs_col = {}
        data = defaultdict(list)
        det_name = flat.getDetector().getName()
        det = camera[det_name]
        response_vs_col[det_name], divisadero_tearing \
            = self.normed_mean_response_vs_col(det, flat.getImage())
        for amp_name, value in divisadero_tearing.items():
            data['det_name'].append(det_name)
            data['amp_name'].append(amp_name)
            data['divisadero_tearing'].append(value)
        return pipeBase.Struct(divisadero_response=np.array(response_vs_col),
                               divisadero_stats=pd.DataFrame(data))

    def normed_mean_response_vs_col(self, det, flat):
        response_vs_col = [[], []]
        for amp in det:
            amp_image = flat.Factory(flat, amp.getBBox())
            amp_data = amp_image.array.copy()
            # We only need to do the flip in y in order to get the
            # rows to consider (via self.config.row_[min,max])
            # correct.  By not flipping in x, we can simply
            # concatenate the mean_by_row results, since the pixels
            # are already in physical order.
            if amp.getRawFlipY():
                amp_data[:] = amp_data[::-1, :]

            # Compute robust mean for each column.
            mean_by_col, _, _ = astropy.stats.sigma_clipped_stats(
                amp_data[self.row_min:self.row_max, :], axis=0)

            # Normalize by median of mean values.
            mean_by_col /= np.median(mean_by_col)

            # Fit a line to mean_by_row and divide by that line.
            nedge = self.nedge
            x = np.arange(mean_by_col.shape[0])
            y = mean_by_col
            cpoly = np.polyfit(x[nedge:-nedge], y[nedge:-nedge], deg=1)
            yfit = cpoly[1] + cpoly[0]*x
            mean_by_col /= yfit

            if amp.getName() in 'C10 C11 C12 C13 C14 C15 C16 C17':
                # "Top" of the CCD (the entire CCD for WFS)
                response_vs_col[0].extend(mean_by_col)
            else:
                # "Bottom" of the CCD.  In the detector object, amps are
                # indexed right to left (C07 -> C06 -> ... -> C00), so
                # flip each amp contribution before extending.
                response_vs_col[1].extend(np.flip(mean_by_col))
        # Convert to numpy arrays.
        response_vs_col[0] = np.array(response_vs_col[0])
        # Flip the bottom-of-CCD array to get left-to-right pixel ordering.
        response_vs_col[1] = np.flip(response_vs_col[1])

        # Loop over leftmost amps of adjacent amp pairs and look for
        # tearing at amp boundaries.
        divisadero_tearing = {}
        # Include a leading edge zero for ghost amp preceding C10.
        max_values = [0]
        # Do the top of CCD first.
        for i in range(7):
            imin = det[i].getBBox().getMaxX() - self.npix
            imax = det[i].getBBox().getMaxX() + 2*self.npix
            max_divisadero \
                = np.nanmax(np.abs(response_vs_col[0][imin:imax] - 1.0))
            if np.isnan(max_divisadero):
                max_divisadero = 0
            max_values.append(max_divisadero)
        # Add trailing edge zero for ghost amp following C17.
        max_values.append(0)
        # Assign values by amp, using maximum of either edge.
        for i, amp_pair in zip(range(8), zip(max_values[:-1], max_values[1:])):
            divisadero_tearing[det[i].getName()] = max(amp_pair)

        # Do bottom of CCD, if it's not a WFS.
        if len(response_vs_col[1]) > 0:
            # The loop starting from C00 counts down from 15 because of
            # the serpentine ordering of amplifiers.
            max_values = [0]
            for i in range(15, 8, -1):
                imin = det[i].getBBox().getMaxX() - self.npix
                imax = det[i].getBBox().getMaxX() + 2*self.npix
                max_divisadero \
                    = np.nanmax(np.abs(response_vs_col[1][imin:imax] - 1.0))
                if np.isnan(max_divisadero):
                    max_divisadero = 0
                max_values.append(max_divisadero)
            max_values.append(0)
            for i, amp_pair in zip(range(15, 7, -1),
                                   zip(max_values[:-1], max_values[1:])):
                divisadero_tearing[det[i].getName()] = max(amp_pair)

        return response_vs_col, divisadero_tearing


class DivisaderoFpPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                       dimensions=("instrument",)):
    divisadero_stats = cT.Input(
        name="divisadero_stats",
        doc=("Divisadero statistics for all CCD amp boundaries "
             "in the focal plane"),
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

    divisadero_tearing = cT.Output(
        name="divisadero_tearing_plot",
        doc="Focal plane mosaic of divisadero tearing results",
        storageClass="Plot",
        dimensions=("instrument",))


class DivisaderoFpPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                                  pipelineConnections=DivisaderoFpPlotsTaskConnections):
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    zmin = pexConfig.Field(doc="Minimum of color bar range.",
                           dtype=float, default=0)
    zmax = pexConfig.Field(doc="Maximum of color bar range.",
                           dtype=float, default=0.05)

class DivisaderoFpPlotsTask(pipeBase.PipelineTask):
    """
    Create a focal plane mosaic of the divisadero tearing results.
    """
    ConfigClass = DivisaderoFpPlotsTaskConfig
    _DefaultName = 'divisaderoFpPlotsTask'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.yfigsize, self.config.xfigsize
        self.z_range = self.config.zmin, self.config.zmax

    def run(self, divisadero_stats, camera):
        # Unpack the divisadero data for plotting.
        amp_data = defaultdict(dict)
        for handle in divisadero_stats:
            det = camera[handle.dataId['detector']]
            det_name = det.getName()
            df = handle.get()
            for _, row in df.iterrows():
                amp_name = row['amp_name']
                amp_data[det_name][amp_name] = row['divisadero_tearing']

        figure = plt.figure(figsize=self.figsize)
        ax = figure.add_subplot(111)
        plot_focal_plane(ax, amp_data, camera=camera, z_range=self.z_range,
                         title='Divisadero max dev', scale_factor='1')

        return pipeBase.Struct(divisadero_tearing=figure)
