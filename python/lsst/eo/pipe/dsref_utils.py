from collections import defaultdict


__all__ = ['RaftOutputRefsMapper']


class RaftOutputRefsMapper:
    def __init___(self, camera):
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
                raft = detector_raft_map[detector]
                ref_map[raft] = ref
        return ref_map
