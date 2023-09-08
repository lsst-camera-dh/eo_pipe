from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import lsst.afw.math as afwMath
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as cT
from .dsref_utils import RaftOutputRefsMapper, get_plot_locations_by_dstype
from .plotting import append_acq_run, nsigma_range


__all__ = ['BiasStabilityTask', 'BiasStabilityPlotsTask']


def get_plot_locations(repo, collections):
    dstypes = ('bias_mean_vs_time_plot', 'bias_stdev_vs_time_plot',
               'bias_rc_mean_vs_time_plot', 'bias_serial_profile_plots',
               'bias_parallel_profile_plots')
    return get_plot_locations_by_dstype(repo, collections, dstypes)


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
    exposures = cT.Input(name="postISRCCD",
                         doc="CCDs with desired level of ISR applied.",
                         storageClass="Exposure",
                         dimensions=("instrument", "exposure", "detector"),
                         multiple=True,
                         deferLoad=True)
    camera = cT.PrerequisiteInput(name="camera",
                                  doc="Camera used in observations",
                                  storageClass="Camera",
                                  isCalibration=True,
                                  dimensions=("instrument",))
    bias_stability_stats = cT.Output(name="bias_stability_stats",
                           doc="bias stability statistics",
                           storageClass="DataFrame",
                           dimensions=("instrument", "detector"))
    bias_serial_profile_plots = cT.Output(
        name="bias_serial_profile_plots",
        doc="profiles of median per-column bias "
        "level as a function of serial pixel",
        storageClass="Plot",
        dimensions=("instrument", "detector"))
    bias_parallel_profile_plots = cT.Output(
        name="bias_parallel_profile_plots",
        doc="profiles of median per-column bias "
        "level as a function of parallel pixel",
        storageClass="Plot",
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
    acq_run = pexConfig.Field(doc="Acquisition run number.",
                              dtype=str, default="")
    xfigsize = pexConfig.Field(doc="Figure size x-direction in inches.",
                               dtype=float, default=16)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=16)
    nsigma_y_axis = pexConfig.Field(doc="Number of sigma to use for y-axis "
                                    "range of profile plots",
                                    dtype=float, default=10)


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
        self.nsigma_y_axis = self.config.nsigma_y_axis
        self.figsize = self.config.xfigsize, self.config.yfigsize

    def run(self, exposures, camera):
        raw_det = camera[exposures[0].dataId['detector']]
        namps = len(raw_det)
        data = defaultdict(list)
        profile_plots = {'serial': plt.figure(figsize=self.figsize),
                         'parallel': plt.figure(figsize=self.figsize)}
        axes = {key: {amp: plot.add_subplot(4, 4, amp)
                      for amp in range(1, namps+1)}
                for key, plot in profile_plots.items()}
        values = {key: {amp: [] for amp in range(1, namps+1)}
                  for key in profile_plots}
        for handle in exposures:
            exp = handle.get()
            image = exp.getImage()
            md = exp.getMetadata()
            det = exp.getDetector()
            det_name = det.getName()
            raft, slot = det_name.split('_')
            # Loop over segments and extract statistics
            for i, amp in enumerate(det, 1):
                amp_name = amp.getName()
                data['run'].append(md['RUNNUM'])
                data['exposure_id'].append(handle.dataId['exposure'])
                data['mjd'].append(md['MJD'])
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
                # Plot the profiles, i.e., the column/row median as a
                # function of x/y-pixel coordinate for this amp.
                imarr = image[bbox].array
                if raw_amp.getRawFlipX():
                    imarr[:] = imarr[:, ::-1]
                y_values = np.median(imarr, axis=0)
                values['serial'][i].extend(y_values)
                axes['serial'][i].plot(range(imarr.shape[1]), y_values)
                if raw_amp.getRawFlipY():
                    imarr[:] = imarr[::-1, :]
                y_values = np.median(imarr, axis=1)
                values['parallel'][i].extend(y_values)
                axes['parallel'][i].plot(range(imarr.shape[0]), y_values)
        title = append_acq_run(self, 'median signal (ADU) vs column', det_name)
        profile_plots['serial'].suptitle(title)
        profile_plots['serial'].tight_layout(rect=(0, 0, 1, 0.95))
        title = append_acq_run(self, 'median signal (ADU) vs row', det_name)
        profile_plots['parallel'].suptitle(title)
        profile_plots['parallel'].tight_layout(rect=(0, 0, 1, 0.95))
        for key, ax in axes.items():
            for i, amp in enumerate(det, 1):
                ax[i].annotate(f"{amp.getName()}", (0.5, 0.95),
                               xycoords='axes fraction', ha='center')
                y_range = nsigma_range(values[key][i],
                                       nsigma=self.nsigma_y_axis)
                if y_range is not None:
                    ax[i].set_ylim(*y_range)
        return pipeBase.Struct(bias_stability_stats=pd.DataFrame(data),
                               bias_serial_profile_plots=profile_plots['serial'],
                               bias_parallel_profile_plots=profile_plots['parallel'])


class BiasStabilityPlotsTaskConnections(pipeBase.PipelineTaskConnections,
                                         dimensions=("instrument",)):
    bias_stability_stats = cT.Input(name="bias_stability_stats",
                                    doc="bias stability statistics",
                                    storageClass="DataFrame",
                                    dimensions=("instrument", "detector"),
                                    multiple=True,
                                    deferLoad=True)
    camera = cT.PrerequisiteInput(name="camera",
                                  doc="Camera used in observations",
                                  storageClass="Camera",
                                  isCalibration=True,
                                  dimensions=("instrument",))
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
                               dtype=float, default=12)
    yfigsize = pexConfig.Field(doc="Figure size y-direction in inches.",
                               dtype=float, default=12)
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
                       'bias stability, mean of region '
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
                    det_name = '_'.join((raft, slot))
                    df = df0.query(f"det_name=='{det_name}'")
                    fig.add_subplot(3, 3, i)
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
                suptitle = self._plot_title_map[plot_type]
                plt.suptitle(append_acq_run(self, suptitle, raft))
                butlerQC.put(fig, ref_map[raft])
                plt.close()
