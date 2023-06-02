import os
from collections import defaultdict
from astropy.io import fits
import matplotlib.pyplot as plt
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT

import scope
import multiscope

from .dsref_utils import RaftOutputRefsMapper, get_plot_locations_by_dstype


__all__ = ['ScanModeTask']


def get_plot_locations(repo, collections):
    dstypes = ('scan_mode_dispersion_plot', 'scan_mode_multiscope_plot')
    return get_plot_locations_by_dstype(repo, collections, dstypes)


def get_raft_arrays(raft_files):
    """Get raftarrrays list for passing to scan mode plotting code."""
    seglist = []
    raftarrays = []
    slots = sorted(list(raft_files.keys()))
    for slot in slots:
        channels = list(range(8)) if slot.startswith('SW') else None
        raftarrays.append(scope.get_scandata_fromfile(raft_files[slot],
                                                      selectchannels=channels))
        seglist.append(slot[1:])
    return raftarrays, seglist


def get_tm_mode(scan_mode_file):
    tm_mode = 'TM_OFF'
    with fits.open(scan_mode_file) as hdus:
        reb_cond = hdus['REB_COND'].header
        if 'AP0_TM' in reb_cond and reb_cond['AP0_TM'] == 1:
            tm_mode = 'TM_ON'
    return tm_mode


class ScanModeTaskConnections(pipeBase.PipelineTaskConnections,
                              dimensions=("instrument", "detector")):
    raw_handles = cT.Input(
        name="raw",
        doc="Raw FITS files from a scan mode frame.",
        storageClass="Exposure",
        dimensions=("instrument", "exposure", "detector"),
        multiple=True,
        deferLoad=True)

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera used in observations",
        storageClass="Camera",
        isCalibration=True,
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibration)

    scan_mode_dispersion_plot = cT.Output(
        name="scan_mode_dispersion_plot",
        doc="CCD-level scan mode dispersion plot.",
        storageClass="Plot",
        dimensions=("instrument", "exposure", "detector"),
        multiple=True)

    scan_mode_multiscope_plot = cT.Output(
        name="scan_mode_multiscope_plot",
        doc="Raft-level scan mode multiscope plot.",
        storageClass="Plot",
        dimensions=("instrument", "exposure", "detector"),
        multiple=True)


class ScanModeTaskConfig(pipeBase.PipelineTaskConfig,
                         pipelineConnections=ScanModeTaskConnections):
    pass


class ScanModeTask(pipeBase.PipelineTask):
    """
    Make scan mode plots.
    """
    ConfigClass = ScanModeTaskConfig
    _DefaultName = "scanModeTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        raw_handles = inputs['raw_handles']
        camera = inputs['camera']
        raft_data = defaultdict(dict)
        acq_run = None
        for ref in raw_handles:
            if acq_run is None:
                acq_run = ref.dataId.records['exposure'].science_program
            raft = ref.dataId.records['detector'].raft
            detector = ref.dataId['detector']
            raft_data[raft][detector] = ref
        self.run(acq_run, raft_data, camera, butlerQC, outputRefs)

    def run(self, acq_run, raft_data, camera, butlerQC, outputRefs):
        # Map the output references for the per-ccd dispersion plots.
        disp_ref_map = {}
        for ref in outputRefs.scan_mode_dispersion_plot:
            det_name = camera[ref.dataId['detector']].getName()
            disp_ref_map[det_name] = ref

        # Map the multiscope output references for each raft.
        raft_output_refs_mapper = RaftOutputRefsMapper(camera)
        multiscope_ref_map = raft_output_refs_mapper.create(
            outputRefs.scan_mode_multiscope_plot)

        # Loop over rafts, retrieve the per-CCD file artifacts for
        # each raft, extract the scan mode data, and make plots.
        for raft, exposure_refs in raft_data.items():
            raft_files = {}
            for detector, ref in exposure_refs.items():
                butler = ref.butler
                det_name = camera[detector].getName()
                slot = det_name.split('_')[1]
                file_name = os.path.basename(butler.getURI(ref).path)
                raft_files[slot] = os.path.join('.', file_name)
            ref.butler.retrieveArtifacts(exposure_refs.values(), "./")
            tm_mode = get_tm_mode(f'./{file_name}')
            # Extract the scan mode data.
            raft_arrays, seg_list = get_raft_arrays(raft_files)
            # Clean up the file artifacts.
            for item in raft_files:
                os.remove(item)
            # Make the dispersion plots for each CCD in the raft.
            for seg, scandata in zip(seg_list, raft_arrays):
                slot = 'S' + seg
                det_name = '_'.join((raft, slot))
                disp_plot_title = f"{det_name}, Run {acq_run}, {tm_mode}"
                Nchan = 8 if slot.startswith('SW') else 16
                scope.plot_scan_dispersion(scandata, title=disp_plot_title,
                                           Nchan=Nchan)
                butlerQC.put(plt.gcf(), disp_ref_map[det_name])
                plt.close()
            # Make the multiscope plot for the current raft.
            title = f"{raft}, Run {acq_run}, {tm_mode}"
            multiscope.plot_raft_allchans(raft_arrays, seg_list, suptitle=title)
            butlerQC.put(plt.gcf(), multiscope_ref_map[raft])
            plt.close()
