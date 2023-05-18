from collections import defaultdict
import pandas as pd


def convert_amp_data_to_df(amp_data_dict, camera=None):
    """
    Convert a triply-nested dict of per-amp measurements, keyed by
    measurement name, det_name, and amp_name to a pandas DataFrame.

    Parameters
    ----------
    amp_data_dict : dict(dict(dict(float)))
       dict containing the per-amp data.
    camera : lsst.afw.cameraGeom.Camera [None]
       Camera object to use to get detector names.

    Returns
    -------
    pandas.DataFrame
    """
    data = defaultdict(list)
    detectors = []
    for column, amp_data in amp_data_dict.items():
        if not data:
            for detector in amp_data:
                if camera is not None:
                    det_name = camera[detector].getName()
                else:
                    det_name = detector
                for amp_name in amp_data[detector]:
                    data['det_name'].append(det_name)
                    data['amp_name'].append(amp_name)
                    detectors.append(detector)
        for det_name, amp_name, detector in zip(data['det_name'],
                                                data['amp_name'],
                                                detectors):
            data[column].append(amp_data[detector][amp_name])
    return pd.DataFrame(data)
