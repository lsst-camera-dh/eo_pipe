from collections import defaultdict
import numpy as np
import pandas as pd
import astropy.stats
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT


__all__ = ['DivisaderoTearingTask']


class DivisaderoTearingTaskConnections(pipeBase.PipelineTaskConnections,
                                       dimensions=("instrument",)):
    flats = cT.Input(
        name="flat",
        doc="Input combined flat to process",
        storageClass="Exposure",
        dimensions=("band", "instrument", "detector", "physical_filter"),
        isCalibration=True,
        multiple=True,
        deferLoad=True)

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
        dimensions=("instrument",))

    divisadero_stats = cT.Output(
        name="divisadero_stats",
        doc=("Divisadero statistics for all CCD amp boundaries "
             "in the focal plane"),
        storageClass="DataFrame",
        dimensions=("instrument",))


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
    _defaultName = 'divisaderoTearingTask'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.row_min = self.config.row_min
        self.row_max = self.config.row_max
        self.nedge = self.config.nedge
        self.npix = self.config.npix

    def run(self, flats, camera):
        response_vs_col = {}
        data = defaultdict(list)
        for handle in flats:
            det = camera[handle.dataId['detector']]
            det_name = det.getName()
            flat = handle.get().getImage()
            response_vs_col[det_name], max_vals \
                = self.normed_mean_response_vs_col(det, flat)
            # For a given amp, the divisadero tearing value is the
            # maximum value of either boundary.  For edge amps, pad
            # with zeros so that the single common boundary is used.
            for amp, devs in zip(det, zip([0] + max_vals, max_vals + [0])):
                data['det_name'].append(det_name)
                data['amp_name'].append(amp.getName())
                data['divisadero_tearing'].append(max(devs))
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
        # Do the top of CCD first:
        divisadero_tearing = []
        for i in range(7):
            imin = det[i].getBBox().getMaxX() - self.npix
            imax = det[i].getBBox().getMaxX() + 2*self.npix
            max_divisadero \
                = np.nanmax(np.abs(response_vs_col[0][imin:imax] - 1.0))
            if np.isnan(max_divisadero):
                max_divisadero = 0
            divisadero_tearing.append(max_divisadero)
        # Do bottom of CCD, if it's not a WFS.
        if len(response_vs_col[1]) > 0:
            # The loop starting from C00 counts down from 15 because of
            # the serpentine ordering of amplifiers.
            for i in range(15, 8, -1):
                imin = det[i].getBBox().getMaxX() - self.npix
                imax = det[i].getBBox().getMaxX() + 2*self.npix
                max_divisadero \
                    = np.nanmax(np.abs(response_vs_col[1][imin:imax] - 1.0))
                if np.isnan(max_divisadero):
                    max_divisadero = 0
                divisadero_tearing.append(max_divisadero)

        return response_vs_col, divisadero_tearing
