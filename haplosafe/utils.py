def complement(base):
    if base == "A":
        return "T"
    elif base == "C":
        return "G"
    elif base == "T":
        return "A"
    elif base == "G":
        return "C"
    else:
        raise ValueError(
            "Base must be 'A', 'C', 'T', or 'G'. Recieved value {} instead.".format(
                base
            )
        )


def reverse_complement(read):
    revcomp = "".join(map(complement, read[::-1].upper()))
    return revcomp


def parse_fasta(filepath):

    with open(filepath) as fasta_file:
        reads = [
            line.strip()
            for i, line in enumerate(fasta_file.readlines())
            if ((i % 4) == 1)
        ]

    return reads


def parse_sam(filepath):
    def _record_check(line):
        return (len(line.split()) > 9) and (len(line.split()[9]) > 100)

    with open(filepath) as sam_file:
        reads = [
            line.split()[9].strip()
            for line in sam_file.readlines()
            if _record_check(line)
        ]

    return reads


def load_reads(filepaths):

    reads = []
    for filepath in filepaths:
        extension = filepath.split(".")[-1]
        if extension == "sam":
            reads.extend(parse_sam(filepath))
        elif extension == "fas":
            reads.extend(parse_fasta(filepath))
        else:
            raise ValueError(
                f"filepath: expected extension to one of '.sam', '.fas', found {extension}"
            )

    # rev_reads = (map(reverse_complement, reads))
    # reads.extend(rev_reads)

    return reads
