from collections import defaultdict
import lsst.daf.butler as daf_butler


__all__ = ['RaftOutputRefsMapper', 'get_plot_locations_by_dstype']


class RaftOutputRefsMapper:
    def __init__(self, camera):
        # Map the representative detector number to raft name.
        det_lists = defaultdict(list)
        for det in camera:
            det_name = det.getName()
            raft, _ = det_name.split('_')
            det_lists[raft].append(det.getId())
        self.detector_raft_map = {min(det_list): raft
                                  for raft, det_list in det_lists.items()}

    def create(self, output_refs):
        # Map the dataset reference for the representative detector
        # to the raft name.
        ref_map = {}
        for ref in output_refs:
            detector = ref.dataId['detector']
            if detector in self.detector_raft_map:
                raft = self.detector_raft_map[detector]
                ref_map[raft] = ref
        return ref_map


def get_plot_locations_by_dstype(repo, collections, dstypes):
    butler = daf_butler.Butler(repo, collections=collections)
    file_paths = defaultdict(list)
    for dstype in dstypes:
        for ref in set(butler.registry.queryDatasets(dstype, findFirst=True)):
            file_paths[dstype].append(butler.getURI(ref).path)
    return dict(file_paths)
