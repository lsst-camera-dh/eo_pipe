from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
from astro_metadata_translator import ObservationInfo
import lsst.afw.math as afwMath
from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from .isr_utils import apply_overscan_correction
from .dsref_utils import RaftOutputRefsMapper
from .plotting import append_acq_run


__all__ = ['BiasStabilityTask', 'BiasStabilityPlotsTask']


def image_stats(image, nsigma=10):
    """Compute clipped mean and stdev of the image."""
    stat_ctrl = afwMath.StatisticsControl(numSigmaClip=nsigma)
    flags = afwMath.MEANCLIP | afwMath.STDEVCLIP
    stats = afwMath.makeStatistics(image, flags=flags, sctrl=stat_ctrl)
    return stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)


def get_readout_corner(image, amp, raw_amp, dxy):
    """
    Return an Image object for the region dxy x dxy pixels including the
    corner where the amp is read out.

    Parameters
    ----------
    image : lsst.afw.image.Image
        Fully assembled raw image from the LSST code.
    amp : lsst.afw.cameraGeom.Amplifier
        Amplifier object for the desired amp using the detector object
        associated with the assembled image.
    raw_amp : lsst.afw.cameraGeom.Amplifier
        Amplifier object for the desired amp for the unassembled image
        so that the flips in x and y required for assembly can be
        ascertained.
    dxy : int
        Size in pixels of the square region to extract that includes
        the imaging pixel at readout corner.

    Returns
    -------
    lsst.afw.image.Image
    """
    bbox = amp.getBBox()
    xmin = bbox.getMaxX() - dxy if raw_amp.getRawFlipX() else bbox.getMinX()
    ymin = bbox.getMaxY() - dxy if raw_amp.getRawFlipY() else bbox.getMinY()
    llc = lsst.geom.Point2I(xmin, ymin)
    extent = lsst.geom.Extent2I(dxy, dxy)
    corner_bbox = lsst.geom.Box2I(llc, extent)
    return image[corner_bbox]


class BiasStabilityTaskConnections(pipeBase.PipelineTaskConnections,
                                   dimensions=("instrument", "detector")):
    raw_frames = cT.Input(name="raw_frames",
                          doc="Raw input frames",
                          storageClass="Exposure",
                          dimensions=("instrument", "exposure", "detector"),
                          multiple=True,
                          deferLoad=True)
    camera = cT.PrerequisiteInput(name="camera",
                                  doc="Camera used in observations",
                                  storageClass="Camera",
                                  isCalibration=True,
                                  dimensions=("instrument",),
                                  lookupFunction=lookupStaticCalibration)
    bias_stability_stats = cT.Output(name="bias_stability_stats",
                           doc="bias stability statistics",
                           storageClass="DataFrame",
                           dimensions=("instrument", "detector"))


class BiasStabilityTaskConfig(pipeBase.PipelineTaskConfig,
                              pipelineConnections=BiasStabilityTaskConnections):
    ncol_skip = pexConfig.Field(doc="Number of initial and trailing columns "
                                "to skip in evaluating serial overscan level",
                                dtype=int, default=4)
    oscan_method = pexConfig.Field(doc="Serial overscan subtraction method",
                                   dtype=str, default="median_per_row")
    readout_corner_size = pexConfig.Field(doc="Size of readout corner region "
                                          "to consider, in pixels.",
                                          dtype=int, default=200)
    nsigma = pexConfig.Field(doc="Number of sigma for [MEAN,STDEV]CLIP",
                             dtype=float, default=10.0)


class BiasStabilityTask(pipeBase.PipelineTask):
    """Task to compute bias stability statistics"""
    ConfigClass = BiasStabilityTaskConfig
    _DefaultName = "biasStabilityTask"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.ncol_skip = self.config.ncol_skip
        self.oscan_method = self.config.oscan_method
        self.readout_size = self.config.readout_corner_size
        self.nsigma = self.config.nsigma

    def run(self, raw_frames, camera):
        raw_det = camera[raw_frames[0].dataId['detector']]
        data = defaultdict(list)
        for handle in raw_frames:
            exp = handle.get()
            image = exp.getImage()
            obs_info = ObservationInfo(exp.getMetadata())
            det = exp.getDetector()
            det_name = det.getName()
            raft, slot = det_name.split('_')
            for amp in det:
                amp_name = amp.getName()
                apply_overscan_correction(exp, amp_name, dx=self.ncol_skip,
                                          method=self.oscan_method)
                data['run'].append(obs_info.science_program)
                data['exposure_id'].append(obs_info.exposure_id)
                data['mjd'].append(obs_info.datetime_begin.value)
                data['det_name'].append(det_name)
                data['amp_name'].append(amp_name)
                bbox = amp.getBBox()
                mean, stdev = image_stats(image[bbox], nsigma=self.nsigma)
                data['mean'].append(mean)
                data['stdev'].append(stdev)
                raw_amp = raw_det[amp_name]
                rc_image = get_readout_corner(image, amp, raw_amp,
                                              self.readout_size)
                rc_mean, rc_stdev = image_stats(rc_image, nsigma=self.nsigma)
                data['rc_mean'].append(rc_mean)
                data['rc_stdev'].append(rc_stdev)
        return pipeBase.Struct(bias_stats=pd.DataFrame(data))


class BiasStabilityPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                         dimensions=("instrument",)):
    bias_stability_stats = cT.Input(name="bias_stability_stats",
                                    doc="bias stability statistics",
                                    storageClass="DataFrame",
                                    dimensions=("instrument", "detector"))
    camera = cT.PrerequisiteInput(name="camera",
                                  doc="Camera used in observations",
                                  storageClass="Camera",
                                  isCalibration=True,
                                  dimensions=("instrument",),
                                  lookupFunction=lookupStaticCalibration)
    bias_mean_vs_time_plot = cT.Output(name="bias_mean_vs_time_plot",
                                       doc="Plot of bias mean vs time",
                                       storageClass="Plot",
                                       dimensions=("instrument", "detector"),
                                       multiple=True)
    bias_stdev_vs_time_plot = cT.Output(name="bias_stdev_vs_time_plot",
                                        doc="Plot of bias stdev vs time",
                                        storageClass="Plot",
                                        dimensions=("instrument", "detector"),
                                        multiple=True)
    bias_rc_mean_vs_time_plot = cT.Output(name="bias_rc_mean_vs_time_plot",
                                          doc="Plot of bias readout corner "
                                          "mean vs time",
                                          storageClass="Plot",
                                          dimensions=("instrument", "detector"),
                                          multiple=True)


class BiasStabilityPlotsTaskConfig(pipeBase.PipelineTaskConfig,
                                   pipelineConnections=BiasStabilityPlotsTaskConnections):
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=9)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=9)
    acq_run = pexConfig.Field(doc="Acquisition run number.",
                              dtype=str, default="")


class BiasStabilityPlotsTask(pipeBase.PipelineTask):
    """
    Task to generate raft-level plots of bias stability statistics as a
    function of time.
    """
    ConfigClass = BiasStabilityPlotsTaskConfig
    _DefaultName = "biasStabilityPlotsTask"

    _plot_column_map = {'bias_mean_vs_time': 'mean',
                        'bias_stdev_vs_time': 'stdev',
                        'bias_rc_mean_vs_time': 'rc_mean'}
    _plot_ylabel_map = {'bias_mean_vs_time': 'mean (ADU)',
                        'bias_stdev_vs_time': 'stdev (ADU)',
                        'bias_rc_mean_vs_time': 'readout corner mean (ADU)'}
    _plot_title_map = {'bias_mean_vs_time': 'bias stability, mean signal',
                       'bias_stdev_vs_time': 'bias stability, stdev',
                       'bias_rc_mean_vs_time':
                       'bias stability, %(self.dxy)sx%(self.dxy)s region '
                       'covering the readout corner'}
    _science_raft_slots = 'S20 S21 S22 S10 S11 S12 S00 S01 S02'.split()
    _corner_raft_slots = 'SG0 SW1 SW0 SG1'.split()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.figsize = self.config.xfigsize, self.config.yfigsize

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        camera = inputs['camera']
        raft_dfs = defaultdict(list)
        for ref in inputs['bias_stability_stats']:
            raft, _ = ref.dataId.records['detector'].full_name.split('_')
            raft_dfs[raft].append(ref.get())
        raft_data = {raft: pd.concat(raft_dfs[raft]) for raft in raft_dfs}
        output_refs = {}
        for plot_type in self._plot_column_map:
            output_refs[plot_type] = eval(f'outputRefs.{plot_type}_plot')

        self.run(raft_data, camera, butlerQC, output_refs)

    def run(self, raft_data, camera, butlerQC, output_refs):
        raft_output_refs_mapper = RaftOutputRefsMapper(camera)
        for plot_type, ref_list in output_refs.items():
            column = self._plot_column_map[plot_type]
            ref_map = raft_output_refs_mapper.create(ref_list)
            for raft, df0 in raft_data.items():
                if raft in 'R00 R04 R40 R44':
                    slots = self._corner_raft_slots
                else:
                    slots = self._science_raft_slots
                mjd0 = int(min(df0['mjd']))
                fig = plt.figure(figsize=self.figsize)
                for i, slot in enumerate(slots, 1):
                    df = df0.query(f"slot=='{slot}'")
                    fig.add_subplot(3, 3, i)
                    det_name = '_'.join((raft, slot))
                    det = camera[det_name]
                    for amp in det:
                        amp_name = amp.getName()
                        my_df = df.query(f"amp_name=='{amp_name}'")
                        plt.scatter(my_df['mjd'] - mjd0, my_df[column], s=2,
                                    label=amp_name)
                    xmin, xmax, _, _ = plt.axis()
                    plt.xlim(xmin, 1.2*(xmax - xmin) + xmin)
                    plt.legend(fontsize='x-small')
                    plt.xlabel(f"MJD - {mjd0}")
                    plt.ylabel(self._plot_ylabel_map[plot_type])
                    plt.title(slot)
                plt.tight_layout(rect=(0, 0, 1, 0.95))
                suptitle = self._plot_title_map[plot_type] % locals()
                plt.suptitle(append_acq_run(self, suptitle, raft))
                butlerQC.put(fig, ref_map[raft])
                plt.close()
