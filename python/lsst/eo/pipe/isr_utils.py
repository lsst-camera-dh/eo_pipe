import numpy as np
import lsst.geom


__all__ = ['apply_minimal_isr']


def get_overscan_region(amp_info, dx):
    """
    Get the bounding box of the serial overscan region omitting dx
    pixels at the leading and trailing edges.
    """
    bbox = amp_info.getRawHorizontalOverscanBBox()
    minY = amp_info.getRawBBox().minY
    height = amp_info.getRawBBox().height
    return lsst.geom.Box2I(lsst.geom.Point2I(bbox.minX + dx, minY),
                           lsst.geom.Extent2I(bbox.width - 2*dx, height))


def apply_overscan_correction(raw, amp, dx=4, method='1d_poly', deg=2):
    """
    Apply serial overcan correction to the specified amp.
    This function will modify the raw pixel data in place.

    It returns the full segment, including pre- and overscan regions,
    after applying the correction to all rows.
    """
    raw_image = raw.getImage()
    det = raw.getDetector()
    amp_info = det[amp]
    full_segment = raw_image.Factory(raw_image, amp_info.getRawBBox())
    overscan = raw_image.Factory(raw_image, get_overscan_region(amp_info, dx))
    if method == 'median':
        full_segment -= np.median(overscan.array)
    elif method == 'mean':
        full_segment -= np.mean(overscan.array)
    elif method in ('median_per_row', '1d_poly'):
        row_vals = np.median(overscan.array, axis=1)
        ny = full_segment.array.shape[0]
        if method == '1d_poly':
            rows = np.arange(ny)
            func = np.poly1d(np.polyfit(rows, row_vals, deg))
            row_vals = func(rows)
        for j in range(ny):
            full_segment.array[j, :] -= row_vals[j]
    return full_segment


def subtract_calib(raw, calib, amp, scale=1):
    """
    Subtract the calibration image for the specifed amp.  The scale option
    is used to multiply the calib image, i.e., if it is a dark and scale is
    the dark time.

    This function modifies the raw pixel data in place.
    """
    if calib is None:
        return
    raw_image = raw.getImage()
    amp_info = raw.getDetector()[amp]
    calib_image = calib.getImage()
    if scale != 1:
        calib_image = calib_image.Factory(calib_image, deep=True)
        calib_image *= scale
    calib_amp_info = calib.getDetector()[amp]
    raw_image[amp_info.getRawDataBBox()] \
        -= calib_image[calib_amp_info.getBBox()]


def apply_minimal_isr(raw, bias, dark, amp, dx=4,
                      oscan_method='median_per_row', deg=2):
    """
    Apply a bare-bones ISR consisting of serial overscan correction,
    combined bias subtraction, combined dark correction.

    This function modifies the raw pixel data in place.

    Return the full amplifier segment.
    """
    full_segment = apply_overscan_correction(raw, amp, dx=dx,
                                             method=oscan_method, deg=deg)
    subtract_calib(raw, bias, amp)
    subtract_calib(raw, dark, amp, scale=raw.getMetadata().get('DARKTIME'))
    return full_segment
