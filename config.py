import os
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "your_email@example.com")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")

REQUEST_DELAY = 0.35 if NCBI_API_KEY else 0.4
MAX_RETRIES = 5
RETRY_BACKOFF_BASE = 2.0

RCRS_ACCESSION = "NC_012920.1"

TARGET_COUNT = 501
BATCH_SIZE = 25
SEARCH_PAGE_SIZE = 200
OVERFETCH_FACTOR = 3

MIN_LENGTH = 16500
MAX_LENGTH = 16700
MAX_N_RATIO = 0.01

INDIA_REGIONS = [
    "Punjab", "Amritsar", "Ludhiana", "Chandigarh",
    "Uttar Pradesh", "UP", "Lucknow", "Kanpur", "Varanasi", "Agra",
    "Delhi", "New Delhi",
    "Maharashtra", "Bombay", "Mumbai", "Pune", "Nagpur",
    "Tamil Nadu", "Chennai", "Madras",
    "Karnataka", "Bangalore", "Bengaluru",
    "West Bengal", "Kolkata", "Calcutta",
    "Gujarat", "Ahmedabad", "Surat",
    "Rajasthan", "Jaipur",
    "Kerala", "Kochi",
    "Andhra Pradesh", "Telangana", "Hyderabad",
    "Bihar", "Patna",
    "Madhya Pradesh", "Bhopal", "Indore",
    "Odisha", "Bhubaneswar",
    "Assam", "Guwahati",
    "Haryana",
    "Jammu", "Kashmir", "Himachal Pradesh",
]
REGION_FILTER = " OR ".join(f'"{r}"[All Fields]' for r in INDIA_REGIONS)


SEARCH_QUERY = (
    '"Homo sapiens"[Organism] AND '
    'mitochondrion[Filter] AND '
    'complete genome[Title] AND '
    'India[Country] AND '
    f'({REGION_FILTER})'
)

WORKDIR = Path(__file__).resolve().parent
CACHE_DIR = WORKDIR / "cache"
OUTPUT_CSV = WORKDIR / "mitovarsitypak_india.csv"
OUTPUT_FASTA = WORKDIR / "mitovarsitypak_india.fasta"
LOG_FILE = WORKDIR / "pipeline.log"

CACHE_DIR.mkdir(exist_ok=True)