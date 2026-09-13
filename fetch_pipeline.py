#!/usr/bin/env python3
"""
fetch_pipeline.py

Fetches human mitochondrial DNA sequences from NCBI Nucleotide (GenBank),
applies strict QC benchmarks, and persists clean sequences/metadata.

Requirements:
    pip install biopython python-dotenv

Environment (.env):
    NCBI_EMAIL=you@example.com
    NCBI_API_KEY=your_ncbi_api_key
"""

import os
import csv
import time
import logging
import random
import socket
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError

from dotenv import load_dotenv
from Bio import Entrez, SeqIO

socket.setdefaulttimeout(60)

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

load_dotenv()

NCBI_EMAIL = os.getenv("NCBI_EMAIL")
NCBI_API_KEY = os.getenv("NCBI_API_KEY")

if not NCBI_EMAIL:
    raise EnvironmentError("NCBI_EMAIL must be set in the environment (.env file).")

Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

DB = "nucleotide"
SEARCH_TERM = (
    'Homo sapiens[Organism] AND mitochondrion[Filter] '
    'AND India[Country] AND 16000:17000[Sequence Length]'
)
RETMAX = 500

# Rate limiting: NCBI allows 10 req/sec with API key, 3 req/sec without.
REQUEST_DELAY = 0.34  # seconds between calls
MAX_RETRIES = 5
BASE_BACKOFF = 1.0  # seconds

# QC thresholds
MIN_LENGTH = 16500
MAX_LENGTH = 16700
MAX_N_FRACTION = 0.01  # 1%

# Post-fetch geographic filter (NCBI's [Country] search term is unreliable
# and can return neighboring/mislabeled records, so we re-verify strictly).
TARGET_COUNTRY = "India"

FASTA_OUTPUT = "mitochondrial_cohort.fasta"
CSV_OUTPUT = "metadata_schema.csv"
LOG_FILE = "pipeline_qc.log"

CSV_FIELDS = ["Accession_ID", "Country", "Region", "Length", "N_Count", "N_Fraction"]

# --------------------------------------------------------------------------
# Logging setup
# --------------------------------------------------------------------------

logger = logging.getLogger("mito_pipeline")
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler(LOG_FILE, mode="a")
file_handler.setFormatter(
    logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
)
logger.addHandler(file_handler)

console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
logger.addHandler(console_handler)


# --------------------------------------------------------------------------
# Resilient Entrez wrapper with exponential backoff
# --------------------------------------------------------------------------

def entrez_call_with_retry(func, **kwargs):
    """
    Executes an Entrez API call with exponential backoff retry logic.
    Retries on HTTPError, URLError, and IncompleteRead.
    """
    attempt = 0
    while attempt < MAX_RETRIES:
        try:
            handle = func(**kwargs)
            time.sleep(REQUEST_DELAY)  # rate-limit buffer after every call
            return handle
        except (HTTPError, URLError, IncompleteRead) as e:
            attempt += 1
            if attempt >= MAX_RETRIES:
                logger.error(f"Max retries exceeded for Entrez call: {e}")
                raise
            backoff = BASE_BACKOFF * (2 ** (attempt - 1)) + random.uniform(0, 0.5)
            logger.warning(
                f"Entrez call failed (attempt {attempt}/{MAX_RETRIES}): {e}. "
                f"Retrying in {backoff:.2f}s."
            )
            time.sleep(backoff)
    raise RuntimeError("Unreachable retry state.")


def search_ncbi():
    """Runs ESearch against NCBI Nucleotide and returns a list of GenBank IDs."""
    logger.info(f"Searching NCBI '{DB}' with term: {SEARCH_TERM}")
    handle = entrez_call_with_retry(
        Entrez.esearch, db=DB, term=SEARCH_TERM, retmax=RETMAX
    )
    record = Entrez.read(handle)
    handle.close()
    id_list = record.get("IdList", [])
    logger.info(f"ESearch returned {len(id_list)} IDs.")
    return id_list


def fetch_record(gb_id):
    """Fetches a single GenBank record (GenBank flat-file format) by ID."""
    handle = entrez_call_with_retry(
        Entrez.efetch, db=DB, id=gb_id, rettype="gb", retmode="text"
    )
    try:
        seq_record = SeqIO.read(handle, "genbank")
    finally:
        handle.close()
    return seq_record


# --------------------------------------------------------------------------
# Parsing helpers
# --------------------------------------------------------------------------

def extract_source_qualifiers(seq_record):
    """
    Extracts Country and Region/geo_loc_name-derived qualifiers from the
    'source' feature. Returns (country, region) strings, defaulting to
    'Unknown' if absent.
    """
    country = "Unknown"
    region = "Unknown"

    for feature in seq_record.features:
        if feature.type == "source":
            qualifiers = feature.qualifiers

            # 'country' qualifier is typically formatted as "Country: Region"
            raw_country = qualifiers.get("country", [None])[0]
            if raw_country:
                if ":" in raw_country:
                    parts = raw_country.split(":", 1)
                    country = parts[0].strip()
                    region = parts[1].strip()
                else:
                    country = raw_country.strip()

            # Fallback / supplement: geo_loc_name (newer GenBank schema)
            if country == "Unknown":
                raw_geo = qualifiers.get("geo_loc_name", [None])[0]
                if raw_geo:
                    if ":" in raw_geo:
                        parts = raw_geo.split(":", 1)
                        country = parts[0].strip()
                        region = parts[1].strip()
                    else:
                        country = raw_geo.strip()

            break  # only one 'source' feature expected

    return country, region


def calculate_n_stats(sequence_str):
    """Returns (n_count, n_fraction) for a given sequence string."""
    seq_upper = sequence_str.upper()
    n_count = seq_upper.count("N")
    length = len(seq_upper)
    n_fraction = (n_count / length) if length > 0 else 1.0
    return n_count, n_fraction


# --------------------------------------------------------------------------
# QC Gate
# --------------------------------------------------------------------------

def apply_qc(accession, sequence_str):
    """
    Applies strict QC benchmarks.
    Returns (passed: bool, reason: str, n_count: int, n_fraction: float).
    """
    length = len(sequence_str)

    if not (MIN_LENGTH <= length <= MAX_LENGTH):
        return False, f"Failed QC: Length {length} outside [{MIN_LENGTH}, {MAX_LENGTH}]", None, None

    n_count, n_fraction = calculate_n_stats(sequence_str)

    if n_fraction > MAX_N_FRACTION:
        return (
            False,
            f"Failed QC: N-content {n_fraction:.4%} exceeds {MAX_N_FRACTION:.0%} threshold "
            f"(N_count={n_count})",
            n_count,
            n_fraction,
        )

    return True, "Passed", n_count, n_fraction


# --------------------------------------------------------------------------
# Output writers
# --------------------------------------------------------------------------

def init_csv_if_needed(path):
    """Writes CSV header if the file does not already exist."""
    file_exists = os.path.isfile(path)
    if not file_exists:
        with open(path, mode="w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
            writer.writeheader()


def append_fasta(path, accession, description, sequence_str):
    with open(path, mode="a") as f:
        f.write(f">{accession} {description}\n")
        # wrap sequence at 70 chars per FASTA convention
        for i in range(0, len(sequence_str), 70):
            f.write(sequence_str[i:i + 70] + "\n")


def append_metadata(path, row_dict):
    with open(path, mode="a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writerow(row_dict)


def get_processed_accessions(csv_path):
    """Reads the CSV and returns a set of already processed Accession IDs."""
    if not os.path.exists(csv_path):
        return set()

    processed = set()
    with open(csv_path, mode="r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            processed.add(row["Accession_ID"])
    return processed


# --------------------------------------------------------------------------
# Main pipeline
# --------------------------------------------------------------------------

def run_pipeline():
    init_csv_if_needed(CSV_OUTPUT)

    # LOAD CHECKPOINT: Get already processed IDs
    processed_ids = get_processed_accessions(CSV_OUTPUT)
    logger.info(f"Found {len(processed_ids)} already processed records. Skipping them.")

    id_list = search_ncbi()
    if not id_list:
        logger.warning("No records found for the given search term. Exiting.")
        return

    passed_count = 0
    failed_count = 0

    for gb_id in id_list:
        accession = gb_id  # fallback label until record is parsed
        try:
            seq_record = fetch_record(gb_id)
            accession = seq_record.id

            # --- APPLY CHECKPOINT: Skip if we already have it ---
            if accession in processed_ids:
                logger.info(f"{accession}: Already processed. Skipping.")
                continue
            # ----------------------------------------------------

            sequence_str = str(seq_record.seq)

            country, region = extract_source_qualifiers(seq_record)

            # --- STRICT GEOGRAPHIC FILTER ---
            # NCBI's [Country] search term can return mislabeled or
            # neighboring-country records (e.g., Nepal, Sri Lanka, Unknown),
            # so re-verify the parsed 'country' qualifier directly.
            if country.lower() != TARGET_COUNTRY.lower():
                reason = f"Failed QC: Country mismatch - {country}"
                logger.info(f"{accession}: {reason}")
                failed_count += 1
                continue
            # ---------------------------------

            passed, reason, n_count, n_fraction = apply_qc(accession, sequence_str)

            if not passed:
                logger.info(f"{accession}: {reason}")
                failed_count += 1
                continue

            # Write clean outputs
            append_fasta(
                FASTA_OUTPUT,
                accession,
                seq_record.description,
                sequence_str,
            )
            append_metadata(
                CSV_OUTPUT,
                {
                    "Accession_ID": accession,
                    "Country": country,
                    "Region": region,
                    "Length": len(sequence_str),
                    "N_Count": n_count,
                    "N_Fraction": f"{n_fraction:.6f}",
                },
            )

            logger.info(f"{accession}: Passed QC and written to output.")
            passed_count += 1

        except Exception as e:
            logger.error(f"{accession}: Failed to process - {type(e).__name__}: {e}")
            failed_count += 1
            continue

    logger.info(
        f"Pipeline complete. Total: {len(id_list)}, Passed: {passed_count}, "
        f"Failed/Dropped: {failed_count}"
    )


if __name__ == "__main__":
    run_pipeline()