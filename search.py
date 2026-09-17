import logging
import config
from ncbi_client import NCBIClient

log = logging.getLogger("search")


def collect_candidate_ids(client: NCBIClient, target_count: int, overfetch_factor=config.OVERFETCH_FACTOR):
    wanted = target_count * overfetch_factor
    ids = []
    retstart = 0
    total = None

    while len(ids) < wanted:
        record = client.esearch(config.SEARCH_QUERY, retstart=retstart, retmax=config.SEARCH_PAGE_SIZE)
        batch_ids = record.get("IdList", [])
        total = int(record.get("Count", 0))
        if not batch_ids:
            break
        ids.extend(batch_ids)
        retstart += len(batch_ids)
        log.info(f"Fetched {len(ids)}/{min(wanted, total)} candidate IDs")
        if retstart >= total:
            break

    if total is not None and total < target_count:
        log.warning(f"NCBI only has {total} matching records, fewer than the target of {target_count}")

    return ids
