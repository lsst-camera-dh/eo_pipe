import os
from collections import defaultdict
import click
import pandas as pd
import lsst.daf.butler as daf_butler


__all__ = ['PtcSelection']


class PtcSelection:
    """
    Class to create a chained colllection containing a down-selected
    set of PTC frames.
    """
    def __init__(self, acq_run, repo, instrument, user=None):
        if user is None:
            self.user = os.environ['USER']
        self.raw_collection = f"{instrument}/raw/all"
        collections = [self.raw_collection]
        self.acq_run = acq_run
        self.butler = daf_butler.Butler(repo, collections=collections,
                                        writeable=True)

    def make_chained_ptc_collection(self, in_collection, stride=5,
                                    flux_keyword="CCOBFLUX", detector=23):
        ptc_chain = f"u/{self.user}/run_{self.acq_run}_{stride:02d}_ptc_chain"
        collection_exists = False
        try:
            self.butler.registry.queryCollections(ptc_chain)
        except daf_butler.registry.MissingCollectionError:
            pass
        else:
            collection_exists = True
            print(f"CHAINED collection {ptc_chain} already exists.")
            if not click.confirm("Overwrite?", default=False):
                print("Aborting.")
                return
        ptc_collection = self.make_tagged_ptc_subset(stride=stride,
                                                     flux_keyword=flux_keyword,
                                                     detector=detector)
        my_chain = set(self.butler.registry.getCollectionChain(in_collection))
        my_chain.remove(self.raw_collection)
        my_chain.add(ptc_collection)
        if collection_exists:
            self.butler.registry.removeCollection(ptc_chain)
        self.butler.registry.registerCollection(
            ptc_chain, type=daf_butler.CollectionType.CHAINED)
        self.butler.registry.setCollectionChain(ptc_chain, my_chain)
        return ptc_chain

    def make_tagged_ptc_subset(self, stride=5, flux_keyword="CCOBFLUX",
                               detector=23):
        """
        Make a TAGGED collection of PTC flat pairs that have been
        down-selected using stride, after ordering by the value in the
        flux_keyword.
        """
        ptc_collection \
            = f'u/{self.user}/run_{self.acq_run}_{stride:02d}_ptc_frames'
        collection_exists = False
        try:
            self.butler.registry.queryCollections(ptc_collection)
        except daf_butler.registry.MissingCollectionError:
            pass
        else:
            collection_exists = True
            print(f"TAGGED collection {ptc_collection} already exists.")
            if not click.confirm("Overwrite?", default=False):
                print("Aborting.")
                return

        df0, ref_dict = self.extract_ptc_md(self.acq_run, detector=detector)

        # Sort target fluxes and down-select
        ccob_fluxes = sorted(list(set(df0[flux_keyword])))[::stride]

        # Loop over target flux values, saving ids for exposure pairs
        exp_ids = set()
        for ccob_flux in ccob_fluxes:
            df = df0.query(f"CCOBFLUX=={ccob_flux}")
            if len(df) == 2:
                exp_ids.update(set(df['exposure']))
        refs = [ref_dict[_] for _ in exp_ids]

        if collection_exists:
            self.butler.registry.removeCollection(ptc_collection)
        self.butler.registry.registerCollection(ptc_collection)
        self.butler.registry.associate(ptc_collection, refs)
        return ptc_collection

    def extract_ptc_md(self, acq_run, detector=None):
        """
        Extract the header data relavent for down-selectng PTC datasets
        created using the CCOB wide beam projector.
        """
        where = (f"exposure.science_program='{acq_run}' "
                 "and exposure.observation_type='flat' "
                 "and exposure.observation_reason='flat'")
        if detector is not None:
            where += f" and detector in ({detector})"
        refs = list(set(self.butler.registry.queryDatasets('raw', where=where)
                        .expanded()))
        data = defaultdict(list)
        ref_dict = {}
        for ref in refs:
            exposure = ref.dataId['exposure']
            ref_dict[exposure] = ref
            md = self.butler.get(ref.makeComponentRef("metadata"))
            exposure_rec = ref.dataId.records['exposure']
            data['exposure'].append(exposure)
            data['CCOBFLUX'].append(md.get('CCOBFLUX'))
            data['day_obs'].append(exposure_rec.day_obs)
            data['seq_num'].append(exposure_rec.seq_num)
        return pd.DataFrame(data), ref_dict
