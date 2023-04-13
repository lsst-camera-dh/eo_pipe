"""
Code for raft-level summary plots of divisadero tearing.
"""

from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np


__all__ = ['divisadero_raft_plots', 'make_divisadero_summary_plot']


def divisadero_raft_plots(butler, acq_run, instrument='LSSTCam'):
    """
    Make all of the divisadero summary plots for each raft in the focal plane.
    """
    camera = butler.get("camera", instrument=instrument)

    raft_data = defaultdict(dict)

    # Get divisadero_response data first, i.e., response vs column
    dsrefs = set(butler.registry.queryDatasets('divisadero_response'))
    div_resp_size = len(dsrefs)
    for dsref in dsrefs:
        det = camera[dsref.dataId['detector']]
        det_name = det.getName()
        raft, slot = det.getName().split('_')
        raft_data[raft][slot] = list(butler.get(dsref).tolist()[det_name])

    # Append divisadero_stats (the max divisadero tearing value
    # at amp boundaries) to each raft/slot item.
    dsrefs = set(butler.registry.queryDatasets('divisadero_stats'))
    assert div_resp_size == len(dsrefs)
    for dsref in dsrefs:
        det = camera[dsref.dataId['detector']]
        det_name = det.getName()
        raft, slot = det.getName().split('_')
        raft_data[raft][slot].append(butler.get(dsref))

    # Make a plot for each raft.
    for raft, data in raft_data.items():
        outfile = f"{acq_run}_{raft}_divisadero_tearing.png"
        make_divisadero_summary_plot(data, title=f"Run {acq_run}, {raft}")
        plt.savefig(outfile)


def make_divisadero_summary_plot(raft_data, title=None, figsize=(20, 20)):
    """
    Make summary plot for a single raft.
    """
    fig = plt.figure(figsize=figsize)
    outer = gridspec.GridSpec(3, 3, wspace=0.3, hspace=0.3)
    nskip_edge = 20

    for i, slot in enumerate(raft_data):
        have_wf_sensor = (len(raft_data[slot]) == 2)
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i],
                                                 wspace=0.1, hspace=0.0)
        for j in range(2):
            ncols = len(raft_data[slot][j])
            xpixval = range(ncols)
            if have_wf_sensor and j == 1:
                continue
            # Use max of max_divisidero_tearing to set the range of plots.
            max_divisadero = (raft_data[slot][-1]['divisadero_tearing']
                              .to_numpy()[j*7:j*7 + 8])
            plot_range = np.nanmax(max_divisadero)

            ax = plt.Subplot(fig, inner[j])
            ax.plot(xpixval[nskip_edge:ncols - nskip_edge],
                    raft_data[slot][j][nskip_edge:ncols - nskip_edge])
            ax.set_xlabel('Col #')
            try:
                ax.set_ylim(1.-plot_range, 1.+plot_range)
            except ValueError as eobj:
                # plot_range is probably inf or NaN because of bad pixel
                # data for this sensor, so just skip this plot.
                print('ValueError:', str(eobj))
                continue
            for k in range(1, 8):
                ax.axvline(x=(ncols//8)*k, color='red', ls='--', alpha=0.2)
            if j == 0 and not have_wf_sensor:
                ax.text(0.025, 0.9, f'{slot}', transform=ax.transAxes)
                ax.text(0.825, 0.05, 'Seg 10-17', transform=ax.transAxes)
            elif j == 1 or have_wf_sensor:
                ax.text(0.825, 0.05, 'Seg 00-07', transform=ax.transAxes)
            fig.add_subplot(ax)

    if title is not None:
        plt.suptitle(title)
    return fig
