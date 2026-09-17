import re
import io
from Bio import SeqIO
import config
from ncbi_client import NCBIClient


def batched(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def tag_region(text: str):
    for region in config.INDIA_REGIONS:
        if re.search(rf"\b{re.escape(region)}\b", text, re.IGNORECASE):
            return region
    return "India (unspecified)"


def extract_pmid(gb_record):
    for ref in gb_record.annotations.get("references", []):
        if getattr(ref, "pubmed_id", None):
            return str(ref.pubmed_id)
    return "NA"


def reference_text(gb_record):
    bits = []
    for ref in gb_record.annotations.get("references", []):
        if getattr(ref, "journal", None):
            bits.append(ref.journal)
    return " ".join(bits)


def fetch_records(client: NCBIClient, ids):
    for batch in batched(ids, config.BATCH_SIZE):
        gb_text = client.efetch_text(batch, rettype="gb", retmode="text")
        for gb_record in SeqIO.parse(io.StringIO(gb_text), "genbank"):
            source_bits = []
            for feature in gb_record.features:
                if feature.type == "source":
                    for values in feature.qualifiers.values():
                        source_bits.extend(values)

            source_text = (
                " ".join(source_bits) + " " +
                gb_record.description + " " +
                reference_text(gb_record)
            )

            yield {
                "accession": gb_record.id,
                "length": len(gb_record.seq),
                "sequence": str(gb_record.seq),
                "region": tag_region(source_text),
                "description": gb_record.description,
                "pubmed_id": extract_pmid(gb_record),
            }
# import re
# import io
# from Bio import SeqIO
# import config
# from ncbi_client import NCBIClient


# def batched(seq, size):
#     for i in range(0, len(seq), size):
#         yield seq[i:i + size]


# def tag_region(text: str):
#     for region in config.INDIA_REGIONS:
#         if re.search(rf"\b{re.escape(region)}\b", text, re.IGNORECASE):
#             return region
#     return "India (unspecified)"


# def fetch_records(client: NCBIClient, ids):
#     for batch in batched(ids, config.BATCH_SIZE):
#         gb_text = client.efetch_text(batch, rettype="gb", retmode="text")
#         for gb_record in SeqIO.parse(io.StringIO(gb_text), "genbank"):
#             source_bits = []
#             for feature in gb_record.features:
#                 if feature.type == "source":
#                     for values in feature.qualifiers.values():
#                         source_bits.extend(values)
#             source_text = " ".join(source_bits) + " " + gb_record.description

#             yield {
#                 "accession": gb_record.id,
#                 "length": len(gb_record.seq),
#                 "sequence": str(gb_record.seq),
#                 "region": tag_region(source_text),
#                 "description": gb_record.description,
#             }
