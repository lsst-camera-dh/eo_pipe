import subprocess
import lsst.daf.butler as daf_butler
from lsst.resources import ResourcePath


__all__ = ['ingest_pd_data', 'copy_ts_data']


def copy_ts_data(src_repo, in_collection, dest_repo, out_collection):
    """
    Copy transmission_sensor data from one repo to another, writing
    out to the latter in the specified run collection.
    """
    src_butler = daf_butler.Butler(src_repo, collections=[in_collection])
    refs = list(set(src_butler.registry.queryDatasets('transmission_sensor',
                                                      findFirst=True)))

    dest_butler = daf_butler.Butler(dest_repo, writeable=True)
    dest_butler.registry.registerCollection(
        out_collection, daf_butler.registry.CollectionType.RUN)

    for ref in refs:
        ts = src_butler.get(ref)
        try:
            dest_butler.put(ts, 'transmission_sensor', ref.dataId,
                            run=out_collection)
        except daf_butler.registry.ConflictingDefinitionError:
            print("ConflictingDefinitionError:", ref.dataId)


def ingest_pd_data(acq_run, instrument='LSSTCam', output_run=None,
                   repo='/repo/ir2'):
    """
    Ingest photodiode .ecsv data for all flat frames in a run.

    Parameters
    ----------
    acq_run : str
        Acquistion run number.
    instrument : str ['LSSTCam']
        Instrument used to take the raw data.
    output_run : str [None]
        Run collection name for the photodiode data. If None, then
        use f"{instrument}/photodiode".
    repo : str ['/repo/ir2']
        Data repository.
    """
    if output_run is None:
        output_run = f"{instrument}/photodiode"

    collections = [f"{instrument}/raw/all", output_run]
    butler = daf_butler.Butler(repo, collections=collections)

    where = (f"exposure.science_program='{acq_run}' "
             "and exposure.observation_type='flat'")
    refs = list(set(butler.registry.queryDatasets("raw", where=where)
                    .expanded()))
    refs = sorted(refs, key=lambda ref: ref.dataId['exposure'])

    # Find unique exposures and an associated reference.
    exposure_refs = {ref.dataId['exposure']: ref for ref in refs}
    print(f"Found {len(exposure_refs)} exposures.", flush=True)

    ingested = []
    for ref in exposure_refs.values():
        raw_path = butler.getURI(ref).path
        pd_path = raw_path[:-len('R22_S11.fits')] + 'photodiode.ecsv'
        pd_uri = f"s3://rubin-sts/{pd_path}"
        resource_path = ResourcePath(pd_uri)
        if not resource_path.exists():
            continue
        command = (f"butler ingest-photodiode --output-run {output_run} "
                   f"{repo} {instrument} {pd_uri}")
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as eobj:
            print(eobj, flush=True)
        else:
            ingested.append(pd_uri)
    print(f"Ingested {len(ingested)} datasets.")
