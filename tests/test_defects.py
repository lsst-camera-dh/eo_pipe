"""
Tests for defectsTask code.
"""
import numpy as np
import lsst.afw.image as afw_image
import lsst.afw.detection as afw_detect
from lsst.obs.lsst import LsstCam
from lsst.ip.isr import Defects
from lsst.eo.pipe.defectsTask import tabulate_defects


def _test_tabulate_defects(detector):
    colthresh = 100
    det = LsstCam.getCamera()[detector]
    amp_bbox = None

    ncol = {}
    npix = {}
    for amp in det:
        if amp_bbox is None:
            amp_bbox = amp.getBBox()  # To get size of imaging region per amp
            nx = amp_bbox.width
            ny = amp_bbox.height
        amp_name = amp.getName()
        ncol[amp_name] = 0
        npix[amp_name] = 0

    value = 1.0
    bbox = det.getBBox()
    image = afw_image.ImageF(bbox.width, bbox.height, 0)

    xmin, xmax = bbox.minX, bbox.maxX
    dx = xmax - xmin

    ymin, ymax = bbox.minY, bbox.maxY
    dy = ymax - ymin

    # Add 10 bad columns to C00 (llc), by filling lower half of the segment
    # over a 10-pixel x-range.
    image.array[ymin:ymin + dy//4, 10:20] += value
    ncol['C00'] += 10

    # Add 20 bad columns to C10 (ulc)
    image.array[3*dy//4: ymax, 10:30] += value
    ncol['C10'] += 20

    # 10x10 region of bad pixels, with half in already-bad columns in C00
    j = ymin + 3*dy//8
    image.array[j:j + 10, 15:25] += value
    npix['C00'] += 50

    # Some isolated bad pixels
    xvals = xmin + np.array((100, 300, 400))
    yvals = ymin + np.array((30, 50, 100, 900, 1600))
    for i in xvals:
        for j in yvals:
            image.array[j, i] += value

    npix['C00'] += len(xvals)*len(yvals)

    # Some bad pixels crossing amp boundaries
    delta_x = 100
    x0 = 3*nx - delta_x//2  # C02 - C03 boundary
    x1 = x0 + delta_x
    image.array[ymin:ymin + 10, x0:x1] += value
    npix['C02'] += 10*(x1 - x0)//2
    npix['C03'] += 10*(x1 - x0)//2

    threshold = afw_detect.createThreshold(value/2.0, 'value')
    footprint_list = afw_detect.FootprintSet(image, threshold).getFootprints()
    defects = Defects.fromFootprintList(footprint_list)

    bad_columns, bad_pixels = tabulate_defects(det, defects,
                                               colthresh=colthresh)

    assert npix == bad_pixels
    assert ncol == bad_columns


def test_tabulate_defects():
    _test_tabulate_defects(27)  # ITL CCD
    _test_tabulate_defects(94)  # e2V CCD


if __name__ == '__main__':
    test_tabulate_defects()
