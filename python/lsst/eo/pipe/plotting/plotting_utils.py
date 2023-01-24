import numpy as np
import lsst.afw.math as afwMath


__all__ =['cmap_range']


def cmap_range(image_array, nsig=5):
    pixel_data = np.array(image_array, dtype=float).flatten()
    stats = afwMath.makeStatistics(pixel_data,
                                   afwMath.STDEVCLIP | afwMath.MEDIAN)
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    vmin = max(min(pixel_data), median - nsig*stdev)
    vmax = min(max(pixel_data), median + nsig*stdev)
    return vmin, vmax
