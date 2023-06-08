import pandas as pd
import scipy.signal
import scipy.stats
import numpy as np
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
import lsst.afw.cameraGeom

__all__ = ['BiasShiftsTask']


def readout_order_arr(arr, amp, verbose=False):
    """Extract the image data from an amp, flipped to match readout order;
       i.e. such that the pixel at (0,0) is the chronological start of readout,
       and (n-1, m-1) is the chronological end.
    Parameters
    ----------
    image : `lsst.afw.image.Image` or `lsst.afw.image.MaskedImage`
        Image containing the amplifier of interest.
    amp : `lsst.afw.cameraGeom.Amplifier`
        Amplifier on image to extract.
    Returns
    -------
    output : `lsst.afw.image.Image`
        Image of the amplifier in the desired configuration.
    """
    flipDic = {
        lsst.afw.cameraGeom.ReadoutCorner.LL: (False, False),
        lsst.afw.cameraGeom.ReadoutCorner.LR: (True, False),
        lsst.afw.cameraGeom.ReadoutCorner.UL: (False, True),
        lsst.afw.cameraGeom.ReadoutCorner.UR: (True, True)
    }
    if verbose:
        print(f"Readout corner: {amp.getReadoutCorner()} - Flip? (X,Y): " +
              str(flipDic[amp.getReadoutCorner()]))

    output = arr
    flipCol, flipRow = flipDic[amp.getReadoutCorner()]
    if flipCol:
        output = output[:, ::-1]
    if flipRow:
        output = output[::-1]
    return output


def get_oscan_rows(exp, amp, nskip=3):
    oscan_arr = readout_order_arr(
        exp[amp.getRawSerialOverscanBBox()].image.array, amp)

    return np.mean(oscan_arr[:, nskip:], axis=1)


def find_mask_runs(mask):
    return np.flatnonzero(np.diff(np.r_[np.int8(0),
                                  mask.view(np.int8),
                                  np.int8(0)])).reshape(-1, 2)


def shift_kernel(window=30):
    kernel = np.concatenate([np.arange(window), np.arange(-window+1, 0)])
    kernel = kernel/np.sum(kernel[:window])
    return kernel


def scan_for_shifts(data, window=30, noise_filter=30,
                    threshold=3, skip_rows=30):
    """
    Looks for shifts in baseline in noisy data using a convolution.

    Returns
    -------
    An array of detected shifts with rows:
    [peak of shift, start of shift area, end of shift area]

    """
    local_noise = np.std(butter_highpass_filter(data, 1/noise_filter, 1))
    shift_conv = np.convolve(data, shift_kernel(window), mode='valid')
    shift_conv = np.concatenate(
        [np.zeros(window-1), shift_conv, np.zeros(window)])
    shift_like = np.abs(shift_conv)/local_noise
    shift_mask = shift_like > threshold
    shift_mask[:skip_rows] = False
    shift_regions = find_mask_runs(shift_mask)
    shift_peaks = []
    for region in shift_regions:
        region_peak = np.argmax(shift_like[region[0]:region[1]]) + region[0]
        if satisfies_flatness(region_peak, shift_conv[region_peak], data):
            shift_peaks.append(
                [shift_conv[region_peak], region_peak, region[0], region[1]])
    if len(shift_peaks) == 0:
        return np.array([[np.NaN]*4]), local_noise
    return np.asarray(shift_peaks),  local_noise


def satisfies_flatness(shiftrow, shiftmag, oscans, window=30, verbose=False):

    prerange = np.arange(shiftrow-window, shiftrow)
    postrange = np.arange(shiftrow, shiftrow+window)

    preslope, preintercept, pre_r_value, p_value, std_err = \
        scipy.stats.linregress(prerange, oscans[prerange])
    postslope, postintercept, post_r_value, p_value, std_err = \
        scipy.stats.linregress(postrange, oscans[postrange])

    if verbose:
        print(f'Pre-threshold: {preslope*2*len(prerange)}')
        print(f'Post-threshold: {postslope*2*len(postrange)}')
        print(f'Shift value: {shiftmag}')

    pretrending = (preslope*2*len(prerange) < shiftmag) \
        if (shiftmag > 0) else (preslope*2*len(prerange) > shiftmag)
    posttrending = (postslope*2*len(postrange) < shiftmag) \
        if (shiftmag > 0) else (postslope*2*len(postrange) > shiftmag)

    return (pretrending and posttrending)


def find_shifts_in_exposure(exp, amp=None, **options):

    if amp is None:
        det = exp.getDetector()
    else:
        det = [amp]

    shifts = {}

    for thisAmp in det:
        oscan_rows = get_oscan_rows(exp, thisAmp)
        shifts[thisAmp.getName()] = scan_for_shifts(oscan_rows, **options), \
            np.mean(oscan_rows), \
            np.mean(exp[thisAmp.getRawDataBBox()].image.array)

    return shifts


def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = scipy.signal.butter(
        order, normal_cutoff, btype='high', analog=False)
    return b, a


def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = scipy.signal.filtfilt(b, a, data)
    return y


class BiasShiftsTaskConnections(pipeBase.PipelineTaskConnections,
                                dimensions=("instrument", "detector")):
    raws = cT.Input(name="raw",
                    doc="Raw input frames",
                    storageClass="Exposure",
                    dimensions=("instrument", "exposure", "detector"),
                    multiple=True,
                    deferLoad=True)

    bias_shifts = cT.Output(name="bias_shifts",
                            doc="Sudden shifts in overscan biases",
                            storageClass="DataFrame",
                            dimensions=("instrument", "detector"))

    bias_shifts_stats = cT.Output(name="bias_shifts_stats",
                                  doc="Per-amp statistics on bias shifts",
                                  storageClass="DataFrame",
                                  dimensions=("instrument", "detector"))


class BiasShiftsTaskConfig(pipeBase.PipelineTaskConfig,
                           pipelineConnections=BiasShiftsTaskConnections):
    pass


class BiasShiftsTask(pipeBase.PipelineTask):
    ConfigClass = BiasShiftsTaskConfig
    _DefaultName = "biasShiftsTask"

    headers = ['Bias Shift (ADU)', 'Shift center (row)', 'Shift start (row)',
               'Shift end (row)', 'HF row noise (ADU)', 'Overscan mean (ADU)',
               'Imaging mean (ADU)', 'raft', 'sensor', 'amp', 'obs_id']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, raws):
        bias_shifts = []
        self.log.info("Found %d raw files", len(raws))

        for i, handle in enumerate(raws):
            exp = handle.get()
            md = exp.getMetadata()
            det = exp.getDetector()
            det_name = det.getName()
            raft, slot = det_name.split('_')

            self.log.info("processing %d: %s", i, handle.dataId['exposure'])

            new_shifts_dic = find_shifts_in_exposure(exp)
            for ampname, ((shifts, hfNoise),
                          oscanMean, imageMean) in new_shifts_dic.items():
                for shift in shifts:
                    bias_shifts.append((shift[0], shift[1], shift[2], shift[3],
                                        hfNoise, oscanMean, imageMean,
                                        raft, slot, ampname, md['OBSID']))

        bias_shifts_df = pd.DataFrame(bias_shifts,
                                      columns=self.headers)
        bias_shifts_df = bias_shifts_df.sort_values(
            by=['obs_id', 'raft', 'sensor', 'amp'])

        shifts_stats = []
        for (raft, sensor, amp), amp_data in \
                bias_shifts_df.groupby(['raft', 'sensor', 'amp']):
            amp_stats = {'raft': raft, 'sensor': sensor, 'amp': amp}
            amp_shifts = amp_data[amp_data['Bias Shift (ADU)'].notna()]

            if len(amp_shifts) == 0:
                amp_stats['shift_freq'] = 0.
                amp_stats['shift_count'] = 0.
                amp_stats['max_shift'] = np.nan
                amp_stats['median_shift'] = np.nan
                amp_stats['min_row'] = np.nan
                amp_stats['max_row'] = np.nan
                amp_stats['median_row'] = np.nan
                amp_stats['std_row'] = np.nan
            else:
                n_obs = amp_data['obs_id'].nunique()
                n_obs_shifts = amp_shifts['obs_id'].nunique()
                amp_stats['shift_freq'] = n_obs_shifts/n_obs

                amp_stats['shift_count'] = len(amp_shifts)
                amp_stats['max_shift'] = np.max(amp_shifts['Bias Shift (ADU)'])
                amp_stats['median_shift'] = np.median(
                    amp_shifts['Bias Shift (ADU)'])
                amp_stats['min_row'] = np.min(amp_shifts['Shift center (row)'])
                amp_stats['max_row'] = np.max(amp_shifts['Shift center (row)'])
                amp_stats['median_row'] = np.median(
                    amp_shifts['Shift center (row)'])
                amp_stats['std_row'] = np.std(amp_shifts['Shift center (row)'])

            shifts_stats.append(amp_stats)
        shifts_stats_df = pd.DataFrame(shifts_stats)
        return pipeBase.Struct(bias_shifts=bias_shifts_df,
                               bias_shifts_stats=shifts_stats_df)
