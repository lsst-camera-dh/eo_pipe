import subprocess
import lsst.daf.butler as daf_butler
from lsst.resources import ResourcePath


__all__ = ['ingest_pd_data', 'copy_repo_data']


def copy_repo_data(src_repo, in_collection, dest_repo, out_collection,
                   dstype='transmission_sensor'):
    """
    Copy transmission_sensor data from one repo to another, writing
    out to the latter in the specified run collection.
    """
    dstype = 'transmission_sensor'
    src_butler = daf_butler.Butler(src_repo, collections=[in_collection])
    refs = list(set(src_butler.registry.queryDatasets(dstype, findFirst=True)))
    if not refs:
        raise RuntimeError(f"Dataset type {dstype} not found in {src_repo} "
                           f"for collection {in_collection}")

    dest_butler = daf_butler.Butler(dest_repo, writeable=True)
    dest_butler.registry.registerCollection(
        out_collection, daf_butler.registry.CollectionType.RUN)

    if not dest_butler.registry.queryDatasetTypes(dstype):
        dest_butler.registry.registerDatasetType(refs[0].datasetType)

    for ref in refs:
        ts = src_butler.get(ref)
        try:
            dest_butler.put(ts, 'transmission_sensor', ref.dataId,
                            run=out_collection)
        except daf_butler.registry.ConflictingDefinitionError:
            print("ConflictingDefinitionError:", ref.dataId)


def ingest_pd_data(acq_run, instrument='LSSTCam', output_run=None,
                   repo='/repo/embargo', s3_bucket="rubin-summit"):
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
    repo : str ['/repo/embargo']
        Data repository.
    s3_bucket : str ['rubin-summit']
        S3 bucket name to use in the data product URI, f"s3://{s3_bucket}/"
    """
    if output_run is None:
        output_run = f"{instrument}/photodiode"

    collections = [f"{instrument}/raw/all", output_run]
    butler = daf_butler.Butler(repo, collections=collections)

    # Find any already-ingested datasets.
    where = f"exposure.science_program='{acq_run}'"
    pd_refs = set(butler.registry.queryDatasets("photodiode", where=where))
    pd_exps = {_.dataId["exposure"] for _ in pd_refs}

    where = (f"exposure.science_program='{acq_run}' "
             "and exposure.observation_type='flat'")
    refs = list(set(butler.registry.queryDatasets("raw", where=where)
                    .expanded()))
    refs = sorted(refs, key=lambda ref: ref.dataId['exposure'])

    # Find unique exposures and an associated reference.
    exposure_refs = {ref.dataId['exposure']: ref for ref in refs
                     if ref.dataId['exposure'] not in pd_exps}
    print(f"Found {len(exposure_refs)} exposures.", flush=True)

    ingested = []
    for ref in exposure_refs.values():
        raw_path = butler.getURI(ref).path
        pd_path = raw_path[:-len('R22_S11.fits')] + 'photodiode.ecsv'
        pd_uri = f"s3://{s3_bucket}/{pd_path}"
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
