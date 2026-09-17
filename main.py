import io
import logging
import traceback
from Bio import SeqIO

import config
from logging_setup import setup_logging
from ncbi_client import NCBIClient
from search import collect_candidate_ids
from fetch import fetch_records
from qc import passes_qc
from align import align_to_rcrs
from export import ResultWriter, load_checkpoint, sequence_hash

setup_logging()
log = logging.getLogger("main")


def fetch_rcrs(client: NCBIClient):
    fasta_text = client.efetch_text([config.RCRS_ACCESSION], rettype="fasta", retmode="text")
    record = next(SeqIO.parse(io.StringIO(fasta_text), "fasta"))
    return str(record.seq)


def get_pubmed_title(client: NCBIClient, pmid, cache):
    if pmid == "NA":
        return "NA"
    if pmid in cache:
        return cache[pmid]
    try:
        summary = client.esummary([pmid], db="pubmed")
        title = summary[0].get("Title", "NA") if summary else "NA"
    except Exception as e:
        log.warning(f"Could not fetch PubMed title for PMID {pmid}: {e}")
        title = "NA"
    cache[pmid] = title
    return title


def main():
    seen_accessions, seen_hashes, accepted_count = load_checkpoint()
    if accepted_count:
        log.info(f"Resuming from checkpoint: {accepted_count} sequences already saved")

    remaining = config.TARGET_COUNT - accepted_count
    if remaining <= 0:
        log.info(f"Checkpoint already has {accepted_count}/{config.TARGET_COUNT} sequences, nothing to do")
        return

    client = NCBIClient()
    writer = ResultWriter(resume=True)
    duplicate_count = 0
    region_counts = {}
    pubmed_title_cache = {}

    try:
        log.info("Fetching rCRS reference sequence")
        rcrs_seq = fetch_rcrs(client)

        log.info(f"Searching NCBI for {remaining} more candidate India mtDNA genomes")
        ids = collect_candidate_ids(client, remaining)
        log.info(f"Collected {len(ids)} candidate accessions")

        if not ids:
            log.error("No candidate accessions returned by NCBI search — check SEARCH_QUERY in config.py")
            return

        for record in fetch_records(client, ids):
            if accepted_count >= config.TARGET_COUNT:
                break

            if record["accession"] in seen_accessions:
                continue

            if not passes_qc(record):
                continue

            seq_hash = sequence_hash(record["sequence"])
            if seq_hash in seen_hashes:
                duplicate_count += 1
                log.info(f"{record['accession']} skipped as duplicate sequence")
                continue
            seen_hashes.add(seq_hash)
            seen_accessions.add(record["accession"])

            record["pubmed_title"] = get_pubmed_title(client, record["pubmed_id"], pubmed_title_cache)

            aligned = align_to_rcrs(rcrs_seq, record["sequence"], record["accession"])
            record["aligned_sequence"] = aligned
            record["raw_sequence"] = record["sequence"]
            writer.write(record)

            accepted_count += 1
            region_counts[record["region"]] = region_counts.get(record["region"], 0) + 1
            log.info(f"[{accepted_count}/{config.TARGET_COUNT}] {record['accession']} accepted "
                      f"({record['region']}, {record['length']} bp, PMID {record['pubmed_id']})")

    except Exception:
        log.error("Pipeline crashed, see traceback below. Partial results are saved and resumable.")
        log.error(traceback.format_exc())

    finally:
        writer.close()
        log.info(f"Wrote {accepted_count} total sequences to {config.OUTPUT_CSV} and {config.OUTPUT_FASTA}")
        log.info(f"Skipped {duplicate_count} duplicate sequences this run")
        log.info("Region breakdown (this run):")
        for region, count in sorted(region_counts.items(), key=lambda x: -x[1]):
            log.info(f"  {region}: {count}")


if __name__ == "__main__":
    main()




