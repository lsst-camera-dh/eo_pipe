from collections import defaultdict
import pandas as pd
import lsst.afw.math as afw_math
import lsst.daf.butler as daf_butler

def extract_bias_and_read_noise(exp, edge_buffer=2):
    det = exp.getDetector()
    data = defaultdict(list)
    for amp in det:
        bbox = amp.getRawSerialOverscanBBox()
        bbox.grow(-edge_buffer)
        overscan = exp.getImage()[bbox]
        flags = afw_math.MEDIAN | afw_math.STDEVCLIP
        stats = afw_math.makeStatistics(overscan, flags)
        data['det_name'].append(det.getName())
        data['amp'].append(amp.getName())
        data['median'].append(stats.getValue(afw_math.MEDIAN))
        data['stdevclip'].append(stats.getValue(afw_math.STDEVCLIP))
    return pd.DataFrame(data=data)


repo = '/sdf/home/j/jchiang/u/cp_pipe_runs/R22_S22_testing/test_repo'
collections = ['LSSTCam/raw/all']

instrument = 'LSSTCam'
acq_run = '13162'
detector = 98
observation_type = 'bias'
observation_reason = 'bias'

where = (f"instrument='{instrument}' and exposure.science_program='{acq_run}' "
         f"and exposure.observation_type='{observation_type}' "
         f"and exposure.observation_reason='{observation_reason}'")

butler = daf_butler.Butler(repo, collections=collections)

dsrefs = list(set(butler.registry.queryDatasets('raw', where=where,
                                                detector=detector)))

dfs = [extract_bias_and_read_noise(butler.get(dsref)) for dsref in dsrefs]
df = pd.concat(dfs)
