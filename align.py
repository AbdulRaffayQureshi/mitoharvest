import io
import shutil
import subprocess
import tempfile
import logging
from Bio import SeqIO, Align

log = logging.getLogger("align")

_mafft_path = shutil.which("mafft")

_aligner = Align.PairwiseAligner()
_aligner.mode = "global"
_aligner.match_score = 2
_aligner.mismatch_score = -1
_aligner.open_gap_score = -5
_aligner.extend_gap_score = -0.5


def _align_with_mafft(rcrs_seq, target_seq, target_id):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as ref_f:
        ref_f.write(f">rCRS\n{rcrs_seq}\n")
        ref_path = ref_f.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tgt_f:
        safe_id = target_id.replace(" ", "_")
        tgt_f.write(f">{safe_id}\n{target_seq}\n")
        tgt_path = tgt_f.name

    result = subprocess.run(
        ["mafft", "--quiet", "--keeplength", "--add", tgt_path, ref_path],
        capture_output=True, text=True, check=True,
    )
    records = list(SeqIO.parse(io.StringIO(result.stdout), "fasta"))
    return str(records[-1].seq)


def _align_with_biopython(rcrs_seq, target_seq):
    alignment = _aligner.align(rcrs_seq, target_seq)[0]
    return str(alignment[1])


def align_to_rcrs(rcrs_seq, target_seq, target_id):
    if _mafft_path:
        try:
            return _align_with_mafft(rcrs_seq, target_seq, target_id)
        except Exception as e:
            log.warning(f"mafft failed for {target_id} ({e}), falling back to pairwise aligner")
    return _align_with_biopython(rcrs_seq, target_seq)
