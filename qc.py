import config


def passes_qc(record: dict) -> bool:
    seq = record["sequence"].upper()
    length = len(seq)

    if not (config.MIN_LENGTH <= length <= config.MAX_LENGTH):
        return False

    n_count = seq.count("N")
    if length == 0 or n_count / length > config.MAX_N_RATIO:
        return False

    valid_bases = set("ACGTN")
    if any(base not in valid_bases for base in set(seq)):
        return False

    return True
